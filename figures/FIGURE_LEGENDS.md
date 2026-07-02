# Figure Legends

**Figure 1. Adjusted association of aneurysm morphometry with aSAH severity.**
Forest plot of proportional-odds ordinal regression odds ratios (per 1 SD
increase in each standardized morphometric index) for the three admission
severity scales (modified Fisher, Hunt & Hess, WFNS). Each estimate is from a
separate single-predictor model adjusted for age, sex, smoking, hypertension,
diabetes, statin use, and aneurysm circulation. Top row, complete-case analysis;
bottom row, multiple imputation (MICE, m = 20). Points are colored by
significance: red, false-discovery-rate (FDR) q < 0.05; blue, nominal p < 0.05;
grey, non-significant. OR > 1 indicates higher severity. Size ratio and aspect
ratio are the most consistent correlates; the size-related indices become
FDR-significant across scales after imputation recovers the ~40% of patients lost
to complete-case analysis.

**Figure 2. Morphology separates into a size axis and a shape axis.**
(A) Principal-component loading plot of the seven morphometric indices. PC1 (the
size axis) is dominated by maximal diameter, size ratio, and volume; PC2 (the
shape axis) by neck width, aspect ratio, and height/width ratio. (B) Scree plot
of variance explained (PC1 48%, PC2 22%). (C) Adjusted ordinal odds ratios for
each component predicting severity: PC1 (size) is associated with WFNS grade,
whereas PC2 (shape) is not, indicating that severity tracks the size dimension of
morphology rather than shape.

**Figure 3. Robustness of the morphology–severity association.**
(A) E-values for the primary indices (multiple-imputation estimates): an
unmeasured confounder would need to be associated with both morphology and
severity by the plotted odds ratio to explain away the point estimate (orange) or
the confidence-interval limit (blue). (B) Proportion of patients presenting with
high-grade Hunt & Hess (3–5) across size-ratio tertiles, showing a threshold
pattern concentrated in the highest tertile (error bars, 95% CI). (C) Missing-
not-at-random sensitivity analysis: the size-ratio → Hunt & Hess odds ratio is
essentially unchanged when imputed grades of patients with missing data are
shifted upward by up to 1.5 severity categories.

**Figure 4. Analytic workflow.**
Schematic of the analysis pipeline: single-center aSAH cohort → seven admission-
imaging morphometric indices and four severity outcomes → multiple imputation for
~22% missingness (with complete-case sensitivity) → primary proportional-odds
ordinal regression (per-SD, covariate-adjusted), with PCA, ridge (ranking only),
and dichotomized logistic models as secondary analyses → Benjamini–Hochberg FDR
control → sensitivity suite (proportional-odds Brant test, MI-then-delete, MNAR
delta-shift, E-values, bootstrap-validated AUC).
