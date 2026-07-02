# Statistical Methods and Model Equations

This document gives the formal specification of every model in the analysis
pipeline, the covariate set, and the estimators used. It is the mathematical
companion to the code walkthrough in `CODE_WALKTHROUGH.md`.

## Notation

For patient $i = 1, \dots, N$:

- $Y_i$ = a clinical severity outcome. Three are ordinal
  (Hunt & Hess $\in\{1..5\}$, WFNS $\in\{1..5\}$, modified Fisher $\in\{1..4\}$);
  GCS $\in\{3..15\}$ is treated as continuous.
- $x_i$ = a single aneurysm morphometric index, standardized to mean 0 / SD 1
  (so every coefficient is "per 1 SD increase").
- $\mathbf{z}_i$ = the vector of prespecified covariates being adjusted for.

## What we adjust for

Every adjusted model conditions on the same covariate block $\mathbf{z}_i$:

$$
\mathbf{z}_i = (\text{age}_i,\ \text{sex}_i,\ \text{smoking}_i,\ \text{HTN}_i,\ \text{DM}_i,\ \text{statin}_i,\ \text{circulation}_i)
$$

- **age** — continuous (years)
- **sex** — Male / Female
- **smoking** — Ever / Never (collapsed from the raw Never/Former/Current to
  avoid empty design-matrix cells and rank deficiency)
- **HTN** — hypertension, Yes / No
- **DM** — diabetes mellitus, Yes / No
- **statin** — statin treatment, Yes / No
- **circulation** — anterior / posterior aneurysm location

Each morphometric index is entered **one at a time** (single-predictor-adjusted
models), not all seven simultaneously. This mirrors the manuscript's design and
avoids the severe collinearity among the size-related indices (see PCA below).

---

## 1. Proportional-odds ordinal logistic regression (primary model)

The three grading scales are ordinal, so they are modeled with the cumulative
(proportional-odds) logit model rather than linear regression. For an outcome
with categories $1, \dots, K$ and cut $k \in \{1, \dots, K-1\}$:

$$
\operatorname{logit}\bigl(P(Y_i \le k)\bigr)
= \log\frac{P(Y_i \le k)}{P(Y_i > k)}
= \theta_k - \bigl(\beta\, x_i + \boldsymbol{\gamma}^\top \mathbf{z}_i\bigr)
$$

- $\theta_1 < \theta_2 < \dots < \theta_{K-1}$ are the ordered intercepts (cutpoints).
- $\beta$ is the log-odds ratio for the morphometric index; a **single** $\beta$
  is shared across all cutpoints — this is the *proportional-odds assumption*.
- $\exp(\beta)$ = odds ratio for being in a **higher** severity category per 1 SD
  increase in the index.

Estimation: maximum likelihood (`MASS::polr`). Wald test for $\beta$:
$z = \hat\beta / \operatorname{SE}(\hat\beta)$, $p = 2\,\Phi(-|z|)$.

### Proportional-odds assumption check

The shared-$\beta$ assumption is tested with the **Brant test**, which refits the
model as $K-1$ separate binary logits and tests whether the slopes are equal:

$$
H_0:\ \beta^{(1)} = \beta^{(2)} = \dots = \beta^{(K-1)}
$$

A small omnibus $p$ means proportional odds is violated. In this cohort it holds
for mFS ($p=0.70$) and WFNS ($p=0.46$) but is violated for HH ($p<0.001$); HH is
therefore additionally reported with a partial proportional-odds model and the
dichotomized logistic model.

---

## 2. Multiple imputation (primary missing-data strategy)

Roughly 22% of grading data and 22% of HTN are missing; complete-case analysis
discards ~40% of patients. Under a **missing-at-random (MAR)** assumption we use
multiple imputation by chained equations (MICE), $m = 20$ imputations.

For each incomplete variable $V_j$ a conditional model is fit on all other
variables and used to draw plausible values:

$$
V_j^{(\text{mis})} \sim P\bigl(V_j \mid V_{-j}\bigr), \qquad j = 1, \dots, p
$$

(predictive mean matching for numeric variables, logistic/polytomous models for
categorical). The target model is fit on each completed dataset
$d^{(1)}, \dots, d^{(m)}$, giving estimates $\hat\beta^{(t)}$ with variances
$U^{(t)}$. Results are combined by **Rubin's rules**:

$$
\bar\beta = \frac{1}{m}\sum_{t=1}^{m}\hat\beta^{(t)}, \qquad
\bar U = \frac{1}{m}\sum_{t=1}^{m} U^{(t)}, \qquad
B = \frac{1}{m-1}\sum_{t=1}^{m}\bigl(\hat\beta^{(t)} - \bar\beta\bigr)^2
$$

$$
T = \bar U + \Bigl(1 + \tfrac{1}{m}\Bigr) B, \qquad
\operatorname{SE}(\bar\beta) = \sqrt{T}
$$

$\bar U$ is the within-imputation variance, $B$ the between-imputation variance,
and $(1+1/m)B$ the penalty for finite $m$. The complete-case model is reported
as the conservative sensitivity analysis.

---

## 3. Logistic regression (dichotomized outcomes)

Severity is dichotomized (HH 3–5, mFS 3–4, WFNS 4–5 = "high"):

$$
\operatorname{logit}\bigl(P(Y_i = 1)\bigr) = \alpha + \beta\, x_i + \boldsymbol{\gamma}^\top \mathbf{z}_i,
\qquad \text{OR} = \exp(\beta)
$$

---

## 4. Principal component analysis (composite morphology index)

The seven indices are standardized and decomposed:

$$
\mathbf{X}_{\text{std}} = \mathbf{U}\mathbf{D}\mathbf{V}^\top,
\qquad \text{PC}_j = \mathbf{X}_{\text{std}}\, \mathbf{v}_j
$$

where $\mathbf{v}_j$ is the $j$-th loading vector and the variance explained by
component $j$ is $d_j^2 / \sum_k d_k^2$. In this cohort:

- **PC1** (48% variance) loads on max diameter, size ratio, volume → a **size axis**.
- **PC2** (22% variance) loads on neck, aspect ratio, H/W ratio → a **shape axis**.

Each component is then used as the predictor $x_i$ in models (1) and (3). Testing
PC2/PC3 (not just PC1) separates the size effect from the shape effect.

---

## 5. Penalized ridge regression (collinearity-robust estimate)

Ridge estimates all indices jointly under an $L_2$ penalty:

$$
\hat{\boldsymbol\beta}^{\text{ridge}}
= \arg\min_{\boldsymbol\beta}\ \Bigl\{ \lVert \mathbf{y} - \mathbf{X}\boldsymbol\beta \rVert_2^2
+ \lambda \lVert \boldsymbol\beta \rVert_2^2 \Bigr\}
$$

with $\lambda$ chosen by 10-fold cross-validation. **Note on inference:** ridge
shrinkage biases the coefficients, so naive Wald $p$-values are not valid. We
therefore report ridge coefficients for *stability/ranking only* and obtain
$p$-values from bootstrap resampling, not from the closed-form covariance. (This
corrects the manuscript's Table 2, which reported invalid ridge $p$-values and
mislabeled the method as "LASSO"; LASSO uses an $L_1$ penalty,
$\lambda\lVert\boldsymbol\beta\rVert_1$, which performs variable selection — a
different method.)

---

## 6. Multiplicity control

Across the primary family of tests we report the Benjamini–Hochberg false
discovery rate. For sorted $p$-values $p_{(1)} \le \dots \le p_{(M)}$:

$$
p^{\text{FDR}}_{(k)} = \min_{j \ge k}\ \min\!\left( \frac{M}{j}\, p_{(j)},\ 1 \right)
$$

Holm and Bonferroni family-wise corrections are reported as stricter sensitivity
checks.

---

## 7. Parsimonious severity index + internal validation

A small multivariable logistic model (size ratio + aspect ratio + bleb) predicts
dichotomized HH. Discrimination is the AUC / c-statistic. Optimism is corrected
with the **Harrell bootstrap** ($B = 500$):

$$
\text{optimism} = \frac{1}{B}\sum_{b=1}^{B}
\Bigl( \text{AUC}(\text{model}_b,\ \text{data}_b) - \text{AUC}(\text{model}_b,\ \text{data}_{\text{orig}}) \Bigr)
$$

$$
\text{AUC}_{\text{corrected}} = \text{AUC}_{\text{apparent}} - \text{optimism}
$$
