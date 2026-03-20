# C160 Final Project Plan: Obesity Paradox — Mechanisms & Variation Analysis (Option A2)

## 1. Project Overview

**Research Question:** Why do obese ICU patients appear to have *lower* hospital mortality than non-obese patients? Decompose the observed mortality difference into direct, indirect, and spurious causal pathways using the Total Variation (TV) decomposition framework from Lecture 14.

**Core Method:** `fairness_cookbook()` from the `faircause` R package (Plečko & Bareinboim, 2024), which implements the TV decomposition theorem:

> TV(x₀, x₁) = x-DE(y | x₀) − x-IE(y | x₀) − x-SE(y)

where DE = direct effect, IE = indirect effect, SE = spurious (confounding) effect.

---

## 2. The Diamond DAG

```
         Z (Age, Sex)
        ↗       ↘
  X (Obese) ——→ Y (Death)
        ↘       ↗
      W (Mediators)
```

| Role | Variable(s) | Justification |
|------|-------------|---------------|
| **X** (Treatment) | `obese` — binarized BMI ≥ 30 | The exposure of interest |
| **Z** (Confounders) | `age`, `sex` | Common causes of both BMI status and mortality. Older patients have higher mortality AND different BMI distributions. Sex affects both obesity prevalence and baseline ICU mortality risk. Neither is *caused* by BMI. |
| **W** (Mediators) | `diag`, `charlson`, `lact_24`, `ast_24`, `pafi_24` | Variables on the causal pathway BMI → W → Mortality. Obesity *causes* certain comorbidities (Charlson), affects illness severity (lactate, PaO2/FiO2), and influences admission type (medical vs. surgical). These then independently affect mortality. |
| **Y** (Outcome) | `death` — hospital mortality flag | Binary indicator of in-hospital death |

**Key DAG justification points for your write-up:**
- Age → Z: Age causes both higher obesity risk and higher mortality; it is not caused by obesity
- Sex → Z: Sex is a biological confounder, not downstream of BMI
- Charlson → W: Obesity causes comorbidities like diabetes, heart failure, which then increase mortality
- Lactate → W: Obese patients may present with different perfusion profiles; lactate reflects acute severity downstream of BMI
- PaO2/FiO2 → W: Obesity affects respiratory mechanics, leading to different oxygenation patterns
- Diagnosis type → W: BMI influences what conditions patients develop, affecting their admission category

---

## 3. Dataset Inventory & Table-by-Table Plan

You confirmed you have all files except `chartevents`. Below is every table needed, what you extract from it, and potential issues.

### 3.1 Files You Have

| File | Size | What You Extract | DAG Role |
|------|------|-----------------|----------|
| `admissions.csv` | ~15 MB | `hospital_expire_flag` → `death`; `hadm_id` for linking | Y (outcome) |
| `icustays.csv` | ~3 MB | `stay_id`, `subject_id`, `hadm_id`, `intime`; restrict to first ICU stay per patient | Cohort definition |
| `omr.csv` | ~15 MB | BMI values → binarize into `obese` (BMI ≥ 30) | X (treatment) |
| `labevents.csv` | ~1.5 GB | Lactate (50813), AST (50878), PaO2 (50821), FiO2 (50816) within 24h of ICU admission | W (mediators) |
| `d_icd_diagnoses.csv` | ~300 KB | **Dictionary only** — maps ICD codes to descriptions. Not the actual diagnosis records. | Reference |
| `d_labitems.csv` | ~40 KB | **Dictionary only** — maps itemid to lab names. Useful for verification but not analysis. | Reference |

### 3.2 Files You Still Need to Download

| File | Source Path | Size | Why You Need It |
|------|------------|------|----------------|
| **`patients.csv.gz`** | `hosp/patients.csv.gz` | ~3 MB | `anchor_age` → `age` and `gender` → `sex` (your Z confounders). **Without this, the entire analysis fails.** |
| **`services.csv.gz`** | `hosp/services.csv.gz` | ~2 MB | `curr_service` → categorized into `diag` (MED/SURG/TRAUM/etc.). This is a W mediator. |
| **`diagnoses_icd.csv.gz`** | `hosp/diagnoses_icd.csv.gz` | ~30 MB | The actual ICD diagnosis codes per admission. Required for computing the Charlson comorbidity index (W mediator). **This is different from `d_icd_diagnoses.csv` which is just the dictionary.** |

**Critical distinction:** `d_icd_diagnoses.csv` ≠ `diagnoses_icd.csv`. The former is a lookup table mapping ICD codes to descriptive names. The latter contains the actual diagnosis codes assigned to each hospital admission, which you need for computing Charlson.

These three files are small (total ~35 MB) and essential. Download them before proceeding.

### 3.3 File You Are Skipping

| File | Why Skip | Impact |
|------|----------|--------|
| `chartevents.csv.gz` | ~3.3 GB compressed, too large | Cannot compute SOFA score. Remove `acu_24` from W. Your mediator set is still strong without it (5 variables). Mention in limitations. |

---

## 4. Code Architecture — What to Fix

Your R code is well-structured. Here are the modifications needed, organized by priority.

### 4.1 Critical Fixes (will cause errors otherwise)

**Fix 1: File paths and extensions.** Your code reads from `~/Downloads/mimiciv/hosp/patients.csv.gz` but your files appear to be flat CSVs. Adjust `DATA_DIR` and file paths to match your actual directory structure:

```r
# Option A: If your files are flat in one folder
DATA_DIR <- "~/Downloads/mimiciv"
patients <- fread(file.path(DATA_DIR, "patients.csv"))  # not .csv.gz, no subdirectory

# Option B: If you keep the MIMIC-IV directory structure
DATA_DIR <- "~/Downloads/mimiciv"
patients <- fread(file.path(DATA_DIR, "hosp", "patients.csv.gz"))
```

Audit every `fread()` call to match your actual file layout.

**Fix 2: Set `COMPUTE_SOFA <- FALSE`.** Since you're skipping `chartevents`, hardcode this to avoid the entire SOFA block:

```r
COMPUTE_SOFA <- FALSE
```

**Fix 3: `diag` as a factor variable.** The `faircause` package's `fairness_cookbook()` uses random forests internally for estimation. Factor variables with many levels can cause issues. You have two options:

- **Option A (recommended):** Convert `diag` to dummy variables before calling `fairness_cookbook()`:
  ```r
  # Create binary indicators for each diagnosis type
  for (dx in unique(dat$diag)) {
    dat[, paste0("diag_", dx) := as.integer(diag == dx)]
  }
  # Drop one reference category (e.g., OTHER) to avoid collinearity
  dat[, diag := NULL]
  W <- c("diag_MED", "diag_SURG", "diag_TRAUM", "diag_NMED", "charlson", 
         "lact_24", "ast_24", "pafi_24")
  ```

- **Option B:** Keep it as a factor but convert to numeric (less interpretable):
  ```r
  dat[, diag := as.numeric(as.factor(diag))]
  ```

Test both and see which runs without error. Option A is cleaner for the write-up.

### 4.2 Analytical Improvements (won't crash, but improve quality)

**Fix 4: BMI timing.** Your current code takes the *most recent* OMR BMI per patient regardless of when it was recorded. Improve by selecting the BMI closest to (but before) the ICU admission:

```r
# Merge with admission time
bmi_data <- merge(bmi_data, 
                  first_stays[, .(subject_id, intime)], 
                  by = "subject_id")
bmi_data[, days_before := as.numeric(difftime(as.Date(intime), chartdate, units = "days"))]
# Keep only BMI recorded before or at ICU admission
bmi_before <- bmi_data[days_before >= 0]
# Take the closest one
bmi_per_patient <- bmi_before[order(subject_id, days_before)][
  , .(bmi = first(bmi_value)), by = subject_id
]
```

If this drops too many patients, fall back to your original approach but document it as a limitation.

**Fix 5: PaO2/FiO2 ratio from same blood gas draw.** Currently you take min(PaO2) and max(FiO2) independently. Better approach: compute the ratio per blood gas draw, then take the worst (lowest) ratio:

```r
# Pair PaO2 and FiO2 by closest timestamp
pao2 <- labs_24[itemid == 50821, .(stay_id, charttime, pao2 = valuenum)]
fio2 <- labs_24[itemid == 50816, .(stay_id, charttime, fio2 = valuenum)]

# Standardize FiO2
fio2[fio2 > 1, fio2 := fio2 / 100]

# For each PaO2 measurement, find the closest FiO2
# (This is the clinically correct approach)
paired <- pao2[fio2, on = .(stay_id, charttime), roll = "nearest"]
paired[, pafi := pao2 / fio2]
pafi_24 <- paired[, .(pafi_24 = min(pafi, na.rm = TRUE)), by = stay_id]
```

If this is too complex, keep your original approach but note the limitation.

**Fix 6: Keep `bmi` for descriptive statistics.** Your code drops the raw BMI column on line 68. Consider keeping it for Table 1 (cohort characteristics), which is expected in any clinical paper. You can drop it just before calling `fairness_cookbook()`.

### 4.3 Presentation Enhancements

**Fix 7: Build a proper heterogeneity heatmap.** Your current code produces separate `autoplot()` calls per subgroup. To match the lecture's style (Lecture 14, slides 23 and 26), create a combined heatmap:

```r
# Collect results into a data.frame
het_results <- data.table()
for (ag in names(het_age)) {
  s <- summary(het_age[[ag]])
  het_results <- rbind(het_results, data.table(
    age_group = ag,
    measure = s$measure,    # adjust to actual column names from summary()
    value = s$value
  ))
}

# Then use ggplot2 with geom_tile() to create the heatmap
```

This is a presentation improvement — do it after the core analysis works.

---

## 5. Step-by-Step Execution Plan

### Phase 1: Data Preparation (do this first)

1. **Download the 3 missing files:** `patients.csv.gz`, `services.csv.gz`, `diagnoses_icd.csv.gz`
2. **Fix file paths** in your code to match your actual directory structure
3. **Set `COMPUTE_SOFA <- FALSE`**
4. **Run Parts 1–10** of your code (data loading through imputation)
5. **Verify cohort size:** You should have roughly 20,000–40,000 patients with valid BMI (OMR coverage is roughly 50–60% of ICU patients)
6. **Print `str(dat)`** and verify all columns are the correct types (numeric for all except possibly `diag`)

### Phase 2: Core Analysis

7. **Handle `diag` factor** (dummy-code it per Fix 3)
8. **Run `fairness_cookbook()`** — this is the core TV decomposition
9. **Inspect output:** `summary(fcb)` should give you TV, DE, IE, SE with confidence intervals
10. **Interpret:** A negative TV means obese patients have *lower* mortality (the "paradox"). The decomposition tells you why:
    - Positive DE → obesity directly *increases* mortality risk (expected)
    - Negative IE → the indirect path through mediators *reduces* the association (obese patients present with different severity profiles)
    - Positive SE → confounders (age) inflate the apparent protective effect (obese ICU patients are younger)

### Phase 3: Heterogeneity Analysis

11. **By age group:** Run `fairness_cookbook()` on each of: 18–49, 50–64, 65–74, 74+
12. **By diagnosis type:** Run on each of: MED, SURG (removing `diag` from W when stratifying)
13. **Optional — by sex:** Run separately for male and female
14. **Compile into heatmap-style table** matching Lecture 14, slides 23/26

### Phase 4: Write-Up

15. **DAG figure** with justification paragraph (why each variable is Z vs. W)
16. **Table 1:** Cohort characteristics stratified by obese vs. non-obese (age, sex, mortality, Charlson, lab values)
17. **Main decomposition plot** (`autoplot(fcb)`)
18. **Heterogeneity heatmaps/tables**
19. **Discussion:** Connect to obesity paradox literature, explain the mechanism
20. **Limitations section:** BMI measurement timing, simplified mediator set (no SOFA), missing data imputation assumptions, binarization threshold (BMI ≥ 30), unmeasured confounders

---

## 6. Expected Results & Interpretation Guide

Based on Lecture 14 slide 20 (MIMIC-IV results from the original Plečko analysis), you should expect something like:

| Component | Expected Sign | Interpretation |
|-----------|---------------|----------------|
| **TV** | Slightly positive or negative (~-1.5% to +0.8%) | The raw mortality difference between obese and non-obese |
| **DE** (Direct) | Positive (~0.7% to 1.4%) | Obesity directly *increases* mortality, holding mediators fixed — this is the "true" causal effect |
| **IE** (Indirect) | Negative (~-1.2% to -2.4%) | Operating through mediators, obesity is associated with *lower* mortality — obese patients may present with lower acuity scores or different comorbidity profiles |
| **SE** (Spurious) | Positive (~0.2% to 0.7%) | Age confounding — obese ICU patients tend to be younger, which creates a spurious association with lower mortality |

The "paradox" is explained by: **the indirect effect through mediators overwhelms the direct effect**, making the total association appear protective. Obese patients admitted to the ICU tend to be younger (SE) and present with different severity profiles (IE), which masks a genuine harmful direct effect of obesity on mortality.

For heterogeneity, Lecture 14 slide 23 shows the **direct effect is driven by medical admissions** and is strongest in the 50–74 age range. Surgical admissions show much weaker direct effects.

---

## 7. Deliverables Checklist

- [ ] Download `patients.csv.gz`, `services.csv.gz`, `diagnoses_icd.csv.gz`
- [ ] Fix file paths and extensions in R code
- [ ] Set `COMPUTE_SOFA <- FALSE`
- [ ] Handle `diag` factor variable (dummy-coding)
- [ ] Run cohort extraction and verify ~20k–40k patients
- [ ] Run main TV decomposition
- [ ] Run heterogeneity by age group
- [ ] Run heterogeneity by diagnosis type
- [ ] Create DAG figure
- [ ] Create Table 1 (cohort characteristics)
- [ ] Create main decomposition plot
- [ ] Create heterogeneity heatmap
- [ ] Write discussion connecting to obesity paradox literature
- [ ] Write limitations section
