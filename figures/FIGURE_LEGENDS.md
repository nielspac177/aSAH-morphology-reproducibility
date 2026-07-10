# Figure Legends

**Figure 1. Analytic workflow.**
Swim-lane schematic of the analysis pipeline: single-center aSAH cohort → seven
admission-imaging morphometric indices and four severity outcomes (Hunt–Hess,
WFNS, modified Fisher [ordinal]; Glasgow Coma Scale [continuous]) → multiple
imputation for ~22% missingness (with complete-case sensitivity) → primary
proportional-odds ordinal regression (per-SD, adjusted for age, sex, smoking,
hypertension, diabetes, statin use, and aneurysm circulation) → Benjamini–Hochberg
FDR control → sensitivity suite (Brant test, MI-then-delete, MNAR delta-shift,
E-values, bootstrap-validated AUC). The primary-model equation and full adjustment
set are shown in the highlighted panel. Ordinal grades are modelled by proportional
odds; admission GCS (lower = worse) by linear regression. CTA, CT angiography; DSA,
digital subtraction angiography.

**Figure 2. Adjusted association of aneurysm morphometry with aSAH severity.**
Table-with-forest of proportional-odds (cumulative) ordinal odds ratios per 1-SD
increase in each standardized morphometric index, each from a separate
single-index model adjusted for age, sex, smoking, hypertension, diabetes, statin
use, and aneurysm circulation, under multiple imputation (MICE, m = 20). Filled
squares are point estimates (marker area proportional to precision) with capped
95% confidence intervals; OR > 1 indicates higher severity. Odds ratios and
Benjamini–Hochberg false-discovery rates (FDR) are listed at right; bold denotes
FDR < 0.05. Complete-case estimates are given in Table 2. Size ratio, aspect
ratio, volume, and — for Hunt–Hess — maximal diameter are the associations
surviving FDR correction. SD, standard deviation.

**Figure 3. Aneurysm morphology separates into a size axis and a shape axis.**
(A) Principal-component loading biplot of the seven morphometric indices; each
index is coloured by the axis on which it loads more strongly (blue, size
axis/PC1: maximal diameter, size ratio, volume, neck/parent ratio; gold, shape
axis/PC2: aspect ratio, height/width ratio, neck width). All seven indices load
positively on PC1, giving a single-sided size axis. (B) Scree plot of variance
explained (PC1 48%, PC2 22%). (C) Adjusted proportional-odds odds ratios per 1 SD
for each component predicting admission severity, adjusted for the same seven
covariates as the primary analysis (age, sex, smoking, hypertension, diabetes,
statin use, and aneurysm circulation). The size axis (PC1) is associated with WFNS
grade (OR 1.35, 95% CI 1.06–1.73; shown bold, red), whereas the shape axis (PC2)
is not. Marker and whisker style as in Figure 2. Full component-wise estimates,
including PC3 and the GCS analysis, are given in Supplementary Table S1.

**Figure 4. Robustness of the morphology–severity associations.**
(A) E-value at the 95% confidence-interval limit for each of the five
FDR-significant associations (multiple-imputation estimates): the minimum
association (risk-ratio scale) that an unmeasured confounder would need with both
morphology and severity to render the result non-significant; greater distance
from 1 indicates greater robustness (for scale, the smoking–lung-cancer E-value is
≈ 10). Point and confidence-interval E-values for all primary indices are
tabulated in Supplementary Table S2. (B) Proportion of patients presenting with
high-grade Hunt–Hess (3–5) across size-ratio tertiles (complete-case for size
ratio, n = 133; error bars, 95% CI). (C) Missing-not-at-random sensitivity
analysis: the size-ratio → Hunt–Hess odds ratio (adjusted for the same seven
covariates as the primary analysis) is essentially unchanged when the imputed
grades of patients with missing data are shifted upward by up to 1.5 severity
categories; the grey band marks the primary (as-imputed) estimate. MNAR, missing
not at random.
