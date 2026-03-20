# Decomposing the Obesity Paradox: A Causal Variation Analysis of BMI and ICU Mortality in MIMIC-IV

## Course
STATS C160/C260 — Causal Inference for Health Data  
UCLA, Winter 2026  
Instructor: Drago Plečko

## Overview
The "obesity paradox" refers to the counterintuitive observation that obese patients in ICU settings appear to have lower hospital mortality than non-obese patients, despite obesity being a well-established risk factor in the general population. This project uses the Total Variation (TV) decomposition framework (Plečko & Bareinboim, 2024) to investigate this phenomenon in MIMIC-IV data.

Rather than estimating a single total causal effect, we decompose the observed mortality difference between obese (BMI ≥ 30) and non-obese ICU patients into three causal pathways:
- **Direct effect (CtfDE):** the effect of BMI on mortality not mediated through illness severity or comorbidities
- **Indirect effect (CtfIE):** the effect operating through mediators such as the Charlson index, lactate, AST, PaO2/FiO2, and admission type
- **Spurious effect (CtfSE):** the portion of the association driven by confounding through age and sex

We further investigate heterogeneity of these effects across age groups and admission types.

## Key Results
- **Cohort:** 19,208 ICU patients (6,721 obese, 12,487 non-obese)
- **Total Variation:** −3.07% (obese patients have lower observed mortality)
- **Decomposition:** TV (−3.07%) = CtfDE (−0.98%) − CtfIE (+0.74%) − CtfSE (+1.35%)
- Age confounding (CtfSE) accounts for a large portion of the apparent protective effect
- The direct protective association is concentrated among older medical patients (ages 65–74) and is largely absent in surgical patients

## Repository Structure

### `code/`
- `C160_Final_Project_code.Rmd` — R Markdown source file containing the full analysis pipeline: data extraction, cohort construction, TV decomposition, and heterogeneity analysis
- `C160_Project_Plan.md` — Project planning document outlining the DAG, data sources, and analysis steps

### `output/`
- `C160_Final_Project_code.html` — Rendered HTML output of the R Markdown file, including all code, console output, and inline figures

### `report/`
- `Stats C160 Final Project Report.docx` — Final project report (8 pages max, 12pt font) covering data description, causal query, DAG justification, methods, results, and limitations

### `appendices/`
Contains all figures referenced in the report, organized as follows:

| File | Appendix | Description |
|------|----------|-------------|
| `main_decomposition.pdf` | Appendix B | Main TV decomposition plot showing TV, CtfDE, CtfIE, and CtfSE for the full cohort |
| `heatmap_direct_effect.pdf` | Appendix C | Heatmap of the direct effect (CtfDE) cross-stratified by age group and admission type (MED vs. SURG) |
| `heatmap_indirect_effect.pdf` | Appendix D | Heatmap of the indirect effect (CtfIE) cross-stratified by age group and admission type |
| `decomp_age_18_49.pdf` | Appendix E | TV decomposition for patients aged 18–49 |
| `decomp_age_50_64.pdf` | Appendix F | TV decomposition for patients aged 50–64 |
| `decomp_age_65_74.pdf` | Appendix G | TV decomposition for patients aged 65–74 |
| `decomp_age_74+.pdf` | Appendix H | TV decomposition for patients aged 74+ |
| `decomp_diag_med.pdf` | Appendix I | TV decomposition for medical (MED) admissions |
| `decomp_diag_surg.pdf` | Appendix J | TV decomposition for surgical (SURG) admissions |
| `decomp_diag_nmed.pdf` | Appendix K | TV decomposition for neuromedical (NMED) admissions |
| `decomp_diag_traum.pdf` | Appendix L | TV decomposition for trauma (TRAUM) admissions |

Note: Appendix A (Table 1: Cohort characteristics) is included in the report itself as a formatted table, not as a separate figure.

## Data
This analysis uses MIMIC-IV v2.2, available at https://physionet.org/content/mimiciv/2.2/.  
Access requires credentialing through PhysioNet. **Raw data files are not included in this repository** per the MIMIC-IV data use agreement.

Tables used: `patients`, `admissions`, `icustays`, `omr`, `services`, `diagnoses_icd`, `labevents`.

## Methods & Tools
- **R packages:** `data.table`, `ggplot2`, `comorbidity`, `faircause`
- **Decomposition:** `fairness_cookbook()` from the `faircause` package (Plečko & Bareinboim, 2024), which implements TV decomposition using random forest-based propensity score estimation
- **Causal model:** Diamond DAG with X = binarized BMI, Z = {age, sex}, W = {diagnosis type dummies, Charlson index, lactate, AST, PaO2/FiO2}, Y = hospital mortality

## References
- Plečko, D. & Bareinboim, E. (2024). Causal Fairness Analysis: A Causal Toolkit for Fair Machine Learning. *Foundations and Trends in Machine Learning*, 17(3), 304–589.
- Johnson, A., Bulgarelli, L., Shen, L., et al. (2023). MIMIC-IV, a freely accessible electronic health record dataset. *Scientific Data*, 10, 1.
