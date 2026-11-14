# adolescent-faers-adr-glp1-sglt2-insulin-disproportionality

Code, derived data and supplementary files for the manuscript

> Telkar A, Telkar A, Javalgikar A, Madanwale N, Ruikar D, Baligar P. *Adverse drug reaction signals for GLP-1 receptor agonists, SGLT2 inhibitors, insulins, and cardiometabolic agents in adolescents with obesity and type 2 diabetes: a four-metric disproportionality analysis of the FDA Adverse Event Reporting System (2021–2025).*

Archive: <https://doi.org/10.5281/zenodo.22767355> (concept DOI, resolves to the latest version).

> **[!IMPORTANT] Corrected pipeline (V19):** This version corrects a methodological error in the prior submission.  
> Ages are no longer imputed; only records with a **confirmed** age of 12–17 years are retained.  
> Drug roles are restricted to **primary (PS) and secondary (SS) suspect** for the primary disproportionality analysis.  
> The corrected adolescent cohort is **122,552 cases** (down from 403,278, which included age-imputed records).

## Contents

| Folder | Contents |
|---|---|
| `code/` | `pipeline_rerun.py` — corrected FAERS pipeline (Tasks A–D); produces all analytic datasets. `train_temporal.py` — chronological ML pipeline, baselines, McNemar tests, bootstrap intervals. `regen_figs.py` — Figs. 1, 3, 4 and the detailed flowchart. |
| `data/` | Corrected derived datasets: `adolescent_cohort.csv` (122,552 adolescent reports, confirmed ages only), `drug_panel_adolescent.csv` (3,196 drug-panel rows, PS/SS/C roles), `drug_event_pairs.csv` (7,799 unique caseid × ingredient × PT triples), `disproportionality_all.csv` (1,791 signal-detection rows), `table1_drug_counts.csv`, `table4_demographics.csv`. MedDRA licence note. |
| `results/` | `pipeline_summary.json` — cohort and signal counts from the corrected pipeline. `acceptance_tests.json` — 5/5 acceptance tests PASSED. `cohort_stats.json` — age distribution. `results_temporal.json` — ML output from `train_temporal.py`. |
| `supplementary/` | Online Resources 1–11 as submitted. `ESM_8.xlsx` — corrected signal detection (1,791 drug–event pairs across all panels and role filters; 309 with N≥3 in broad obesity PS/SS). |
| `figures/` | Manuscript figures (EPS and PNG). |

## Reproduce the corrected disproportionality pipeline

```bash
pip install -r requirements.txt
# Requires raw FAERS ASCII files in data/raw/ (see FAERS download instructions below)
python code/pipeline_rerun.py   # produces data/ and results/ output files (~2–5 min)
```

**Corrected cohort (from `acceptance_tests.json`, all 5/5 PASSED):**

| Check | Value | Status |
|---|---|---|
| Most common age ≤25% of cohort | 16.0 yr at 18.4% | PASS |
| Cases outside [12,18) | 0 | PASS |
| Duplicate (caseid, ingredient, PT) triples | 0 | PASS |
| All 2×2 cell sums = 122,552 | ✓ | PASS |
| Min signal N ≥ 3 | ✓ | PASS |

## Reproduce the ML classification

```bash
python code/train_temporal.py     # writes results/results_temporal.json (about 20 s)
python code/regen_figs.py         # writes Figs. 1, 3, 4 and the detailed flowchart to figures/
```

`train_temporal.py` reads `data/ObesityAll_14_Drugs_Adolescent_analytic.xlsx`; pass `--data` to use another file and `--out` to write elsewhere.

Expected output with XGBoost 3.2.0, seed 42 (2025 test set, N = 2,812):

| Metric | Value |
|---|---|
| XGBoost temporal Macro-F1 | 0.389 |
| Preferred Term lookup baseline Macro-F1 | 0.594 |
| XGBoost without the Preferred Term | 0.222 |
| McNemar, XGBoost vs lookup | b = 134, c = 453, chi-square = 172.27 |

The XGBoost figure depends on the library version (0.424 with XGBoost 3.4.1), so the version is pinned in `requirements.txt`. Bootstrap interval bounds vary slightly with NumPy and scikit-learn versions.

## Key corrected numbers

| Metric | Prior (flawed) | Corrected |
|---|---|---|
| Raw FAERS records | 7,612,804 | **8,598,529** |
| After deduplication | — | **7,349,591** |
| Excluded (no valid age) | — | **3,189,002** |
| Adolescent cohort (12–17 yr confirmed) | 403,278 | **122,552** |
| Age mean ± SD | 14.57 ± 0.63 (impossible) | **14.69 ± 1.67** |
| Broad Obesity signals (PS/SS, N≥3, ROR) | 360 | **231** |
| Broad Diabetes signals (PS/SS, N≥3, ROR) | — | **151** |
| Metformin PS/SS cases | many | **0 (removed)** |

## FAERS raw data

Raw FAERS files are public: <https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html>  
Download quarterly ASCII files for 2021Q1–2025Q1 and place in `data/raw/`.

## Not included

- Raw FAERS files (public, see above)
- MedDRA. It is licensed from the MSSO and must not be redistributed; see `data/MEDDRA_LICENSE_NOTE.txt`. Fig. 2 cannot be regenerated without a MedDRA licence.
- The fitted model object (needed for the SHAP figure).

## Licence

Code: MIT (see `LICENSE`). Derived data and documents: CC BY 4.0.
