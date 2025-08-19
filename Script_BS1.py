#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script B.S1 (clean version)
------------------------------------
Computes:
  (1) per-sample δ18O_env (Table B.1)
  (2) paleoelevation matrix by subset × model × proxy (Table B.2)

Method summary:
- Carbonates: VPDB→VSMOW; calcite–water equilibrium at T = 20 ± 2.5 °C (Kim & O'Neil, 1997).
- Glass: δD_vg → δD_w via α_vg/w = 0.9668 ± 0.0005 (Friedman, 1993); then LMWL δD = a*δ18O + b with a=8.29, b=11.75.
- Curves: empirical + Rayleigh (winter/summer), all shifted −0.6 ‰ (Miocene baseline).
- Summer curve vertically anchored to winter by matching minimum elevation (no assumption of Miocene elevation)
- Model envelopes (low/high) represent ±2σ; elevation jitter uses σ_model = (z_high − z_low)/4.
- Monte Carlo seed fixed for reproducibility.

Expected inputs (same directory as this script):
  Carbonate_Data_Raw.csv
  Glass_Data_Raw.csv
  lapse_curves.csv      # long format with columns: d18O_VSMOW, z_m OR z_km, curve

Outputs (written to ./outputs):
  outputs/d18Oenv_MC_estimates.csv
  outputs/paleoelev_matrix_summary_95_68CI.csv
"""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

# -------------------------- Paths & I/O helpers --------------------------

def base_dir() -> Path:
    try:
        return Path(__file__).resolve().parent
    except NameError:
        return Path.cwd()

BASE = base_dir()
CARB_FILE  = BASE / "Carbonate_Data_Raw.csv"
GLASS_FILE = BASE / "Glass_Data_Raw.csv"
LAPSE_FILE = BASE / "lapse_curves.csv"

OUT_DIR      = BASE / "outputs"
OUT_DIR.mkdir(exist_ok=True)
OUT_ENV_FILE = OUT_DIR / "d18Oenv_MC_estimates.csv"
OUT_MATRIX   = OUT_DIR / "paleoelev_matrix_summary_95_68CI.csv"

# -------------------------- Global settings -----------------------------

RANDOM_SEED = 42

# Monte Carlo
N_DRAWS_PER_SAMPLE = 100_000    # per-sample δ18O_env draws
N_DRAWS_PER_SUBSET = 100_000    # elevation draws from subset means

# Constants & conversions
T_MEAN_C = 20.0
T_SD_C   = 2.5
ALPHA_VG_W_MEAN = 0.9668
ALPHA_VG_W_SD   = 0.0005
LMWL_SLOPE = 8.29
LMWL_INTERCEPT = 11.75
MIOCENE_D180_SHIFT = -0.6  # permil shift applied to ALL curves

# Subsets (age in Ma; inclusive). For 18–19 Ma carbonate split, use a δ18O_env threshold:
CARB_1819_SPLIT_THRESHOLD = -10.0

SUBSETS = [
    ("33–37 Ma",                 "Carbonate", 33.0, 37.0, None),
    ("17.9–19 Ma (all)",         "Both",      17.9, 19.0, None),
    ("17.9–19 Ma (low-δ18O)",    "Carbonate", 17.9, 19.0, (-np.inf, CARB_1819_SPLIT_THRESHOLD)),
    ("17.9–19 Ma (high-δ18O)",   "Carbonate", 17.9, 19.0, (CARB_1819_SPLIT_THRESHOLD,  np.inf)),
    ("16.5–17.8 Ma",             "Both",      16.5, 17.8, None),
    ("14.1–16.4 Ma",             "Both",      14.1, 16.4, None),
]

# Carbonate filters
INCLUDE_CARBONATE_TYPES = {"Nodule", "Lens", "Rhizolith"}  # exclude Spar
CARBONATE_SOURCE_FILTER = None  # e.g., "ThisStudy" to restrict

# -------------------------- Math utilities -----------------------------

def vpdb_to_vsmow(delta_vpdb: np.ndarray) -> np.ndarray:
    """Coplen (1983) δ18O conversion VPDB → VSMOW."""
    return 1.03091 * delta_vpdb + 30.91

def celsius_to_kelvin(t_c: np.ndarray) -> np.ndarray:
    return t_c + 273.15

def alpha_calcite_water_KimONeil1997(T_kelvin: np.ndarray) -> np.ndarray:
    """1000 ln α = 18030/T − 32.42 (T in K)."""
    return np.exp((18030.0 / T_kelvin - 32.42) / 1000.0)

def delta_water_from_calcite(delta_c_vsmow: np.ndarray, alpha_cw: np.ndarray) -> np.ndarray:
    """δw = [ (δc/1000 + 1)/α − 1 ]*1000."""
    return ((delta_c_vsmow / 1000.0 + 1.0) / alpha_cw - 1.0) * 1000.0

def delta_water_from_glassD(deltaD_glass: np.ndarray, alpha_vg_w: np.ndarray) -> np.ndarray:
    """δD_w = [ (δD_g/1000 + 1)/α_vg/w − 1 ]*1000."""
    return ((deltaD_glass / 1000.0 + 1.0) / alpha_vg_w - 1.0) * 1000.0

def d18O_from_dD_LMWL(deltaD_w: np.ndarray, slope: float = LMWL_SLOPE, intercept: float = LMWL_INTERCEPT) -> np.ndarray:
    """δD = a*δ18O + b  ⇒  δ18O = (δD − b)/a."""
    return (deltaD_w - intercept) / slope

def qtiles(x: np.ndarray, qs=(0.50, 0.16, 0.84, 0.025, 0.975)) -> Dict[str, float]:
    q = np.quantile(x, qs)
    return {"median": float(q[0]), "p16": float(q[1]), "p84": float(q[2]), "p2.5": float(q[3]), "p97.5": float(q[4])}

# -------------------------- Curve container & I/O ----------------------

@dataclass
class Curve:
    """Isotope–elevation curve family with 2 sigma envelopes."""
    d18O_main: np.ndarray; z_main_km: np.ndarray
    d18O_low:  np.ndarray; z_low_km:  np.ndarray
    d18O_high: np.ndarray; z_high_km: np.ndarray
    family:    str  # 'empirical' | 'winter' | 'summer'

def _read_lapse_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    need = {"d18O_VSMOW", "z_km", "curve"}
    if not need.issubset(df.columns):
        missing = need - set(df.columns)
        raise ValueError(f"{path} missing required column(s): {sorted(missing)}")
    return df

def load_curve(path: Path, family: str, shift_d18O: float = 0.0) -> Curve:
    """
    Load <family>-{main,low,high} from a long-format lapse_curves.csv.
    Applies a horizontal δ18O shift (permil) to all three envelopes.
    """
    df = _read_lapse_csv(path)

    def pick(which: str) -> pd.DataFrame:
        tag = f"{family}-{which}"
        sub = df.loc[df["curve"].astype(str).eq(tag), ["d18O_VSMOW", "z_km"]].copy()
        if sub.empty:
            raise ValueError(f"No rows for curve='{tag}' in {path}")
        return sub.sort_values("d18O_VSMOW")

    m  = pick("main"); lo = pick("low"); hi = pick("high")

    d18O_m = m["d18O_VSMOW"].to_numpy() + shift_d18O
    d18O_l = lo["d18O_VSMOW"].to_numpy() + shift_d18O
    d18O_h = hi["d18O_VSMOW"].to_numpy() + shift_d18O

    z_m_km = m["z_km"].to_numpy()
    z_l_km = lo["z_km"].to_numpy()
    z_h_km = hi["z_km"].to_numpy()

    return Curve(d18O_m, z_m_km, d18O_l, z_l_km, d18O_h, z_h_km, family)

def interp_z_from_curve(curve: Curve, d18O_env: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Interpolate z(δ18O) on main/low/high (linear extrapolation at ends)."""
    z_main = np.interp(d18O_env, curve.d18O_main, curve.z_main_km)
    z_low  = np.interp(d18O_env, curve.d18O_low,  curve.z_low_km)
    z_high = np.interp(d18O_env, curve.d18O_high, curve.z_high_km)
    return z_main, z_low, z_high

def draw_elevation_with_model_uncertainty(curve: Curve, d18O_draws: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """
    For each δ18O draw, elevation ~ Normal(z_main, σ_model), with σ_model from ±2σ envelopes:
      σ_model = (z_high − z_low)/4.
    """
    z_m, z_l, z_h = interp_z_from_curve(curve, d18O_draws)
    sigma_model = np.maximum(np.abs(z_h - z_l) / 4.0, 1e-9)
    return rng.normal(loc=z_m, scale=sigma_model)

def reanchor_summer_to_winter_minZ(summer: Curve, winter: Curve) -> Curve:
    """
    Vertical anchor: keep δ18O arrays unchanged; shift all summer z by Δz so
    min(z_summer_main) == min(z_winter_main). Apply same Δz to low/high.
    """
    zmin_w = float(np.min(winter.z_main_km))
    zmin_s = float(np.min(summer.z_main_km))
    delta_z = zmin_w - zmin_s  # often negative (shift summer down)
    return Curve(
        summer.d18O_main, summer.z_main_km + delta_z,
        summer.d18O_low,  summer.z_low_km  + delta_z,
        summer.d18O_high, summer.z_high_km + delta_z,
        summer.family
    )

# -------------------------- Data ingestion -----------------------------

def read_carbonates(path: Path) -> pd.DataFrame:
    """
    Required columns (case-insensitive variants handled):
      Sample_ID, Age_Ma, Type, Source (optional), d18O_mean_VPDB, d18O_sd_VPDB
    Defaults: if d18O_sd_VPDB missing, uses 0.10 ‰.
    Filters: keeps only INCLUDE_CARBONATE_TYPES and optional CARBONATE_SOURCE_FILTER.
    """
    df = pd.read_csv(path)
    colmap = {}
    for c in df.columns:
        cl = c.lower()
        if cl in {"sample_id","sample","id"}: colmap[c] = "Sample_ID"
        elif "age" in cl and "ma" in cl:       colmap[c] = "Age_Ma"
        elif cl == "type":                     colmap[c] = "Type"
        elif "source" in cl:                   colmap[c] = "Source"
        elif "d18" in cl and "vpdb" in cl and any(k in cl for k in ("mean","avg","mu")):
            colmap[c] = "d18O_mean_VPDB"
        elif "d18" in cl and "vpdb" in cl and any(k in cl for k in ("sd","std","sigma")):
            colmap[c] = "d18O_sd_VPDB"
    df = df.rename(columns=colmap)

    need = {"Sample_ID","Age_Ma","Type","d18O_mean_VPDB"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"Carbonate file missing columns: {sorted(missing)}")

    if "d18O_sd_VPDB" not in df.columns:
        df["d18O_sd_VPDB"] = 0.10  # ‰ default

    # type/source filters
    df = df[df["Type"].isin(INCLUDE_CARBONATE_TYPES)].copy()
    if CARBONATE_SOURCE_FILTER is not None and "Source" in df.columns:
        df = df[df["Source"].eq(CARBONATE_SOURCE_FILTER)].copy()

    return df.reset_index(drop=True)

def read_glass(path: Path) -> pd.DataFrame:
    """
    Required columns (case-insensitive variants handled):
      Sample_ID, Age_Ma, dD_vg_VSMOW, dD_sd
    Defaults: if dD_sd missing, uses 1.0 ‰.
    """
    df = pd.read_csv(path)
    colmap = {}
    for c in df.columns:
        cl = c.lower()
        if cl in {"sample_id","sample","id"}: colmap[c] = "Sample_ID"
        elif "age" in cl and "ma" in cl:       colmap[c] = "Age_Ma"
        elif ("dd" in cl or "d2h" in cl) and "vsmow" in cl and any(k in cl for k in ("mean","avg","mu","")):
            colmap[c] = "dD_vg_VSMOW"
        elif ("dd" in cl or "d2h" in cl) and any(k in cl for k in ("sd","std","sigma")):
            colmap[c] = "dD_sd"
    df = df.rename(columns=colmap)

    need = {"Sample_ID","Age_Ma","dD_vg_VSMOW"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"Glass file missing columns: {sorted(missing)}")

    if "dD_sd" not in df.columns:
        df["dD_sd"] = 1.0  # ‰ default

    return df.reset_index(drop=True)

# -------------------------- δ18O_env computations ---------------------

def compute_d18Oenv_carbonate(df: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    for _, r in df.iterrows():
        d18O_c_vpdb = rng.normal(r["d18O_mean_VPDB"], r["d18O_sd_VPDB"], size=N_DRAWS_PER_SAMPLE)
        d18O_c_vsmow = vpdb_to_vsmow(d18O_c_vpdb)
        T_draws = rng.normal(T_MEAN_C, T_SD_C, size=N_DRAWS_PER_SAMPLE)
        alpha_cw = alpha_calcite_water_KimONeil1997(celsius_to_kelvin(T_draws))
        d18O_w = delta_water_from_calcite(d18O_c_vsmow, alpha_cw)
        rows.append({
            "Sample_ID": r["Sample_ID"], "Age_Ma": r["Age_Ma"],
            "Proxy": "Carbonate", "Type": r.get("Type", ""),
            "d18O_env_mean": float(d18O_w.mean()),
            "d18O_env_sd":   float(d18O_w.std(ddof=1))
        })
    return pd.DataFrame(rows)

def compute_d18Oenv_glass(df: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    for _, r in df.iterrows():
        dD_g = rng.normal(r["dD_vg_VSMOW"], r["dD_sd"], size=N_DRAWS_PER_SAMPLE)
        alpha_vg_w = rng.normal(ALPHA_VG_W_MEAN, ALPHA_VG_W_SD, size=N_DRAWS_PER_SAMPLE)
        dD_w = delta_water_from_glassD(dD_g, alpha_vg_w)
        d18O_w = d18O_from_dD_LMWL(dD_w, slope=LMWL_SLOPE, intercept=LMWL_INTERCEPT)
        rows.append({
            "Sample_ID": r["Sample_ID"], "Age_Ma": r["Age_Ma"],
            "Proxy": "Glass", "Type": "Tuff",
            "d18O_env_mean": float(d18O_w.mean()),
            "d18O_env_sd":   float(d18O_w.std(ddof=1))
        })
    return pd.DataFrame(rows)

# -------------------------- Subsets & elevation matrix ----------------

def build_subset_index_split_by_proxy(d_env: pd.DataFrame) -> Dict[Tuple[str, str], np.ndarray]:
    """
    Returns masks keyed by (subset_name, proxy_name).
    'Both' subsets are split into ('Carbonate','Glass'); proxy-specific subsets stay as-is.
    Optional δ18O_env ranges are applied after proxy filtering.
    """
    idx: Dict[Tuple[str,str], np.ndarray] = {}
    for name, flag, amin, amax, drange in SUBSETS:
        age_mask = (d_env["Age_Ma"] >= amin) & (d_env["Age_Ma"] <= amax)
        proxies = ("Carbonate","Glass") if flag == "Both" else (flag,)
        for p in proxies:
            mask = age_mask & d_env["Proxy"].eq(p)
            if drange is not None:
                lo, hi = drange
                mask = mask & (d_env["d18O_env_mean"] >= lo) & (d_env["d18O_env_mean"] <= hi)
            idx[(name, p)] = mask.to_numpy()
    return idx

def draw_subset_mean(d_env: pd.DataFrame, mask: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """
    Draws of the subset-mean δ18O_env: for each sample i, draw Normal(μ_i, σ_i) and average across samples.
    """
    sel = d_env.loc[mask, ["d18O_env_mean", "d18O_env_sd"]].to_numpy()
    if sel.size == 0:
        return np.array([])
    mu = sel[:,0]; sd = sel[:,1]
    draws = rng.normal(mu[:,None], sd[:,None], size=(mu.size, N_DRAWS_PER_SUBSET))
    return draws.mean(axis=0)

# -------------------------- Main --------------------------------------

def main():
    rng = np.random.default_rng(RANDOM_SEED)

    # Load samples and compute per-sample δ18O_env
    carb = read_carbonates(CARB_FILE)
    glas = read_glass(GLASS_FILE)

    d_carb = compute_d18Oenv_carbonate(carb, rng)
    d_glass = compute_d18Oenv_glass(glas, rng)
    d_env = pd.concat([d_carb, d_glass], ignore_index=True)
    d_env.to_csv(OUT_ENV_FILE, index=False)

    # Load curves (apply −0.6 ‰ to ALL)
    empirical = load_curve(LAPSE_FILE, family="empirical", shift_d18O=MIOCENE_D180_SHIFT)
    winter    = load_curve(LAPSE_FILE, family="winter",    shift_d18O=MIOCENE_D180_SHIFT)
    summer    = load_curve(LAPSE_FILE, family="summer",    shift_d18O=MIOCENE_D180_SHIFT)
    # Vertical anchor: match summer min(z) to winter min(z)
    summer    = reanchor_summer_to_winter_minZ(summer, winter)

    # Subsets → elevation matrix (split-by-proxy for 'Both')
    idx = build_subset_index_split_by_proxy(d_env)

    records: List[Dict[str, object]] = []
    for (subset_name, proxy_name), mask in idx.items():
        d_means = draw_subset_mean(d_env, mask, rng)
        if d_means.size == 0:
            continue
        # Elevation posteriors by model
        z_emp = draw_elevation_with_model_uncertainty(empirical, d_means, rng)
        z_win = draw_elevation_with_model_uncertainty(winter,    d_means, rng)
        z_sum = draw_elevation_with_model_uncertainty(summer,    d_means, rng)

        for model_name, z in [("Empirical (mixed-season)", z_emp),
                              ("Rayleigh (Winter)",        z_win),
                              ("Rayleigh (Summer)",        z_sum)]:
            qq = qtiles(z)
            records.append({
                "Subset": subset_name, "Proxy": proxy_name, "Model": model_name,
                "N_samples": int(mask.sum()),
                "z_median_km": qq["median"], "z_p16_km": qq["p16"], "z_p84_km": qq["p84"],
                "z_p2.5_km": qq["p2.5"], "z_p97.5_km": qq["p97.5"],
            })

    matrix = pd.DataFrame.from_records(records)

    # Order rows using categoricals—no helper columns needed
    subset_cats = [name for name,_,_,_,_ in SUBSETS]
    proxy_cats  = ["Carbonate","Glass"]
    model_cats  = ["Rayleigh (Winter)", "Empirical (mixed-season)", "Rayleigh (Summer)"]

    matrix = (matrix
              .assign(Subset=pd.Categorical(matrix["Subset"], categories=subset_cats, ordered=True),
                      Proxy =pd.Categorical(matrix["Proxy"],  categories=proxy_cats,  ordered=True),
                      Model =pd.Categorical(matrix["Model"],  categories=model_cats,  ordered=True))
              .sort_values(["Subset","Proxy","Model"])
              .reset_index(drop=True))

    matrix.to_csv(OUT_MATRIX, index=False)
    print(f"Saved:\n  {OUT_ENV_FILE}\n  {OUT_MATRIX}")

if __name__ == "__main__":
    main()
