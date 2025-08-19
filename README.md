# Script B.S1 — δ18O_env & Paleoelevation Matrix

This repository contains **Script_BS1.py**, which reproduces the supplemental tables:
- **Table D.1:** per-sample estimates of environmental water δ¹⁸O (δ18O_env)
- **Table D.2:** paleoelevation matrix by subset × proxy × model (median and 68/95% CIs)

The implementation follows the manuscript Methods:
- **Carbonates:** VPDB→VSMOW; calcite–water equilibrium at **20 ± 2.5 °C** (Kim & O’Neil, 1997).
- **Glass:** δD_vg → δD_w via **α_vg/w = 0.9668 ± 0.0005** (Friedman, 1993); then **LMWL** (δD = 8.29·δ¹⁸O + 11.75).
- **Curves:** empirical + seasonal Rayleigh (winter/summer), all shifted **−0.6‰** (Miocene baseline).
- **Summer anchoring:** vertical shift so **min(z\_summer) = min(z\_winter)** (x unchanged).
- **Model uncertainty:** curve envelopes treated as **±2σ**, so σ\_model = (z\_high − z\_low)/4.
- **Monte Carlo:** reproducible random seed (**42**); 100,000 draws per sample/subset.

---

## File Layout

```

/ (project root)
├─ \_code/
│  └─ Script\_BS1.py
├─ data/
│  ├─ Carbonate\_Data\_Raw\.csv       # added post-embargo
│  ├─ Glass\_Data\_Raw\.csv           # added post-embargo
│  └─ lapse\_curves.csv             # long format: d18O\_VSMOW, z\_km (or z\_m), curve
└─ outputs/                        # created by the script

````

> **Embargo note:** The raw data files (`Carbonate_Data_Raw.csv`, `Glass_Data_Raw.csv`) will be posted here **after the dissertation embargo is lifted**. The script already uses **relative paths**, so it will run unchanged once the files are added.

---

## Running

1. **Environment**
   - Python ≥ 3.9  
2. **Execute**
   ```bash
   python _code/Script_BS1.py
````

3. **Outputs**

   * `outputs/d18Oenv_MC_estimates.csv`
   * `outputs/paleoelev_matrix_summary_95_68CI.csv`

The script automatically creates `outputs/` and writes results there.

---

## Inputs

* `data/lapse_curves.csv` — long-format curves with columns:

  * `d18O_VSMOW` (‰), `z_km` (preferred; if absent, `z_m` is accepted and treated as km), `curve`
  * `curve` values include: `empirical-{main,low,high}`, `winter-{main,low,high}`, `summer-{main,low,high}`

* `data/Carbonate_Data_Raw.csv` — expects (case-insensitive variants handled):

  * `Sample_ID`, `Age_Ma`, `Type`, `Source` (optional), `d18O_mean_VPDB`, `d18O_sd_VPDB`
  * Defaults: if SD missing, **0.10‰** is used; Spar excluded; Nodules/Lenses/Rhizoliths included.

* `data/Glass_Data_Raw.csv` — expects:

  * `Sample_ID`, `Age_Ma`, `dD_vg_VSMOW`, `dD_sd`
  * Default: if SD missing, **1.0‰** is used.

---

## Reproducibility

* Random seed fixed to **42**.
* All paths are **relative to the script’s location**.
* Curve shift and anchoring are hard-coded to match the manuscript.

---

## Citation

Please cite the dissertation (and archived code/data DOI, once released).
A software/data DOI (e.g., Zenodo) will be added here **post-embargo**.

---

## License

* Code: MIT
* Data: CC BY 4.0

```
