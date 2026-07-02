# Code Walkthrough — Line by Line

A teaching companion to the analysis scripts. The goal is that a reader who does
not know R can follow *exactly* what each line does and *why* it is there. We
quote the code verbatim and explain it block by block. Equations referenced as
"(Eq. N)" point to `METHODS_AND_EQUATIONS.md`.

Pipeline order:

```
00_inspect.R   -> discover how the raw workbook is coded (diagnostic, run once)
01_clean.R     -> build the analysis dataset; verify it against manuscript Table 1
02_analysis.R  -> corrected + extended models; writes outputs/analysis_log.txt
```

---

## Part I — `01_clean.R`: from raw workbook to analysis dataset

### Loading the reader

```r
suppressMessages({
  library(readxl)
})
```

`readxl` reads `.xlsx` without Java or external drivers. `suppressMessages()`
hides the "Attaching package" chatter so the script's own output stays clean.

### Finding the data, keeping PHI out of the repo

```r
xlsx <- Sys.getenv("ASAH_XLSX")
if (!nzchar(xlsx)) stop("Set ASAH_XLSX to the raw workbook path.")
```

The raw workbook contains MRN and DOB (identifiable PHI), so its path is **never
hard-coded**. It is passed through the environment variable `ASAH_XLSX`.
`nzchar()` is `TRUE` only for a non-empty string; if the variable is unset the
script stops immediately rather than failing later with a confusing error.

### Reading past the two-row header

```r
raw <- suppressWarnings(suppressMessages(
  read_excel(xlsx, sheet = 1, skip = 1, .name_repair = "minimal")))
raw <- as.data.frame(raw)
```

The workbook's row 1 is a **category band** ("Baseline characteristics",
"Aneurysm characteristics", …) and row 2 is the real header. `skip = 1` drops the
band so row 2 becomes the column names. `.name_repair = "minimal"` tells `readxl`
not to auto-rename blank/duplicate headers (which would otherwise produce
`...2, ...3, …`). We convert the tibble to a plain `data.frame` for base-R
indexing.

### Two small helper functions

```r
num  <- function(x) suppressWarnings(as.numeric(trimws(as.character(x))))
blank_na <- function(x) { x <- trimws(as.character(x)); x[x %in% c("", "NA", "N/A", "na")] <- NA; x }
```

Several numeric columns (Volume, WFNS, parent-artery diameter) are stored as
**text** in Excel. `num()` trims whitespace and coerces to numeric; the
`suppressWarnings` hides the expected "NAs introduced by coercion" for true
blanks. `blank_na()` normalizes the many spellings of "missing" (empty string,
the literal text `NA`, `N/A`) into a real `NA`, because Excel export leaves
missing values as inconsistent strings.

### Pulling columns by position

```r
d <- data.frame(
  id          = seq_len(nrow(raw)),
  sex_raw     = blank_na(raw[[7]]),
  age         = num(raw[[8]]),
  ...
  aspect_ratio= num(raw[[42]]),   # dome height / neck width -- NOT "ASPECTS"
  ...
)
```

Columns are selected by **integer position** (`raw[[7]]`, `raw[[8]]`, …), not by
name, because the header text contains embedded line breaks and trailing spaces
that make name-matching fragile. The `# NOT "ASPECTS"` comment records a
correction: column 42 is the aneurysm **aspect ratio** (dome height ÷ neck
width), a standard morphology index — *not* the Alberta Stroke Program Early CT
Score, which the manuscript incorrectly expanded.

### Harmonizing the covariates

```r
d$sex <- with(d, ifelse(grepl("^f", sex_raw, ignore.case = TRUE) | sex_raw == "1",
                        "Female",
                 ifelse(grepl("^m", sex_raw, ignore.case = TRUE) | sex_raw == "0",
                        "Male", NA)))
```

The workbook mixes two coding schemes for sex — text (`Female`/`Male`) **and**
numeric (`0`/`1`). `grepl("^f", …)` matches anything starting with "f" (Female,
female, F); the `| sex_raw == "1"` also catches the numeric code. Anything else
becomes `NA`. Nested `ifelse()` is R's vectorized if/else.

```r
smk <- tolower(d$smoke_raw)
d$smoking <- NA_character_
d$smoking[smk %in% c("0", "never smoker", "never tobacco user", "unknown if ever smoked", "never")] <- "Never"
d$smoking[smk %in% c("1", "former", "former smoker", "former tobacco user")] <- "Former"
d$smoking[grepl("^2", smk) | smk %in% c("current", "current smoker")] <- "Current"
d$smoking_ever <- ifelse(is.na(d$smoking), NA, ifelse(d$smoking == "Never", "Never", "Ever"))
```

Smoking had **nine** different codings. We start everyone at `NA` and then assign
`Never`/`Former`/`Current` by matching against every observed spelling.
`grepl("^2", smk)` catches all the "2…" variants ("2", "2 tobacco user").
`smoking_ever` is a two-level collapse used in the adjusted models: with only
~171 patients, a three-level smoking factor creates near-empty design cells that
make the ordinal model **rank-deficient**; Ever/Never avoids that.

```r
bin01 <- function(x) { x <- as.character(x); ifelse(x == "1", "Yes", ifelse(x == "0", "No", NA)) }
d$dm     <- bin01(d$dm_raw)
d$htn    <- bin01(d$htn_raw)
d$statin <- bin01(d$statin_raw)
d$circulation <- ifelse(d$circ_raw == "1", "Posterior",
                 ifelse(d$circ_raw == "0", "Anterior", NA))
```

`bin01()` turns the numeric 0/1 flags into readable `No`/`Yes` labels (anything
else → `NA`). Circulation is mapped 0 → anterior, 1 → posterior.

### Defining the outcomes

```r
d$HH_ord   <- factor(d$HH,   levels = 1:5, ordered = TRUE)
d$mFS_ord  <- factor(d$mFS,  levels = 1:4, ordered = TRUE)
d$WFNS_ord <- factor(d$WFNS, levels = 1:5, ordered = TRUE)
```

`ordered = TRUE` is the crucial line for the whole reanalysis: it tells R the
grades are **ordinal**, which is what lets `MASS::polr` fit the proportional-odds
model (Eq. 1) instead of pretending the grades are continuous.

```r
d$HH_hi   <- ifelse(d$HH   >= 3, 1L, 0L)   # HH 3-5 vs 1-2
d$mFS_hi  <- ifelse(d$mFS  >= 3, 1L, 0L)   # mFS 3-4 vs 1-2
d$WFNS_hi <- ifelse(d$WFNS >= 4, 1L, 0L)   # WFNS 4-5 vs 1-3
```

The dichotomized versions (Eq. 3) use the manuscript's thresholds. `1L` is an
integer literal (the `L` suffix), which keeps these as clean 0/1 integers.

### Transforming skewed predictors and standardizing

```r
d$log_volume <- log(d$volume + 1)
morph <- c("neck", "maxdiam", "volume", "log_volume",
           "aspect_ratio", "hw_ratio", "size_ratio", "neck_parent")
for (m in morph) d[[paste0("z_", m)]] <- as.numeric(scale(d[[m]]))
```

Volume is extremely right-skewed (mean 257, SD 687), so a raw coefficient is
uninterpretable and non-normal; `log(volume + 1)` compresses the tail (the `+1`
guards against `log(0)`). The loop then creates a **z-scored** copy of every
index (`scale()` subtracts the mean and divides by the SD). Standardization is
what makes every odds ratio "per 1 SD", so effects across indices on very
different scales become directly comparable.

### Saving a de-identified dataset

```r
dir.create("../outputs", showWarnings = FALSE)
saveRDS(d, "../outputs/analysis_data.rds")
```

`d` never received MRN, DOB, or any date column, so the saved `.rds` is
de-identified by construction. `showWarnings = FALSE` keeps it quiet if the
folder already exists.

### The verification block (why this script can be trusted)

```r
cat(sprintf("Smoking:  ")); print(table(d$smoking, useNA="ifany"))
cat("   [Never 92, Former 22, Current 51, miss 6]\n")
```

The tail of the script prints each cleaned variable next to the manuscript's
published Table 1 counts (in brackets). If a recoding rule were wrong, the counts
would diverge. They match (smoking 92/22/51, HTN 49/84/38, all morphology means
exact), which validates every mapping above. This is a *reproducibility
self-test*, not decoration.

---

## Part II — `02_analysis.R`: the corrected and extended models

### Setup

```r
suppressMessages({
  library(MASS); library(ordinal); library(mice); library(brant); library(pROC)
})
d <- readRDS("../outputs/analysis_data.rds")
COVARS <- c("age", "sex", "smoking_ever", "htn", "dm", "statin", "circulation")
INDICES <- c("z_size_ratio", "z_aspect_ratio", "z_maxdiam", "z_log_volume",
             "z_neck", "z_hw_ratio", "z_neck_parent")
OUTS <- list(mFS = "mFS_ord", HH = "HH_ord", WFNS = "WFNS_ord")
```

`COVARS` is the adjustment set $\mathbf{z}$ from Eq. (Adjust). `INDICES` are the
seven **standardized** predictors. `OUTS` maps a short label to each ordinal
outcome column so the loops can name their output.

### Block 1 — adjusted ordinal regression, every index × every outcome

```r
for (on in names(OUTS)) {
  for (pr in INDICES) {
    f <- as.formula(paste0(OUTS[[on]], " ~ ", pr, " + ", paste(COVARS, collapse = " + ")))
    m <- tryCatch(MASS::polr(f, data = d, Hess = TRUE, method = "logistic"),
                  error = function(e) NULL)
    if (is.null(m)) next
    ct <- coef(summary(m)); if (!pr %in% rownames(ct)) next
    b <- ct[pr, 1]; z <- ct[pr, 3]; p <- 2 * pnorm(abs(z), lower.tail = FALSE)
    res <- rbind(res, data.frame(outcome = on, predictor = sub("^z_", "", pr),
                 N = nrow(model.frame(m)), OR = exp(b), ...))
  }
}
res$p_FDR <- p.adjust(res$p, "BH")
```

The two nested loops build one **single-predictor-adjusted** model per
(outcome, index) pair. `paste0(...)` assembles a formula string like
`mFS_ord ~ z_size_ratio + age + sex + smoking_ever + htn + dm + statin + circulation`,
and `as.formula()` turns the string into a model formula. `MASS::polr` fits the
proportional-odds model (Eq. 1); `Hess = TRUE` computes the Hessian so standard
errors are available. `tryCatch(..., error = function(e) NULL)` keeps the loop
alive if one model fails to converge. From the coefficient table we read the
estimate `b`, form the Wald `z`, and get the two-sided `p`. `exp(b)` converts the
log-odds to an **odds ratio**. `p.adjust(res$p, "BH")` applies Benjamini–Hochberg
FDR (Eq. 6) across the whole family.

### Block 2 — the same models under multiple imputation

```r
imp <- mice(dm_imp, m = 20, method = meth, printFlag = FALSE, seed = 2024)
```

`mice()` creates 20 completed datasets (Eq. 2). `seed = 2024` fixes the random
draws so the run is reproducible. `printFlag = FALSE` silences the per-iteration
log.

```r
pool_polr <- function(f, pr, imp) {
  ests <- c(); vars <- c()
  for (i in 1:imp$m) {
    di <- complete(imp, i)
    m <- tryCatch(MASS::polr(f, data = di, Hess = TRUE), error = function(e) NULL)
    if (is.null(m)) next
    ct <- coef(summary(m))
    ests <- c(ests, ct[pr, 1]); vars <- c(vars, ct[pr, 2]^2)
  }
  qbar <- mean(ests); ubar <- mean(vars); b <- var(ests)
  T <- ubar + (1 + 1/length(ests)) * b
  se <- sqrt(T); p <- 2 * pnorm(abs(qbar / se), lower.tail = FALSE)
  c(OR = exp(qbar), lo = exp(qbar - 1.96*se), hi = exp(qbar + 1.96*se), p = p)
}
```

This is **Rubin's rules by hand** (Eq. 2). We refit the target model on each
completed dataset `complete(imp, i)`, collect the 20 estimates (`ests`) and their
variances (`vars`). `qbar` is the pooled estimate; `ubar` is the within-imputation
variance $\bar U$; `b` is the between-imputation variance $B$; and
`T = ubar + (1 + 1/m) * b` is the total variance. We do it manually because
`mice::pool()`'s generic does not reliably dispatch on `polr` fits. The pooled
Wald `p` follows as before.

### Block 3 — is the proportional-odds assumption safe?

```r
bt <- brant::brant(m)
cat(sprintf("  %s ~ size_ratio: Omnibus Brant p = %.3f\n", on, bt["Omnibus", "probability"]))
```

`brant()` runs the Brant test (Eq. 1, assumption check). A small omnibus $p$
means the single shared slope is not appropriate. Here HH fails ($p<0.001$), so
HH results are qualified and cross-checked against the dichotomized logistic
model.

### Block 4 — separating the size axis from the shape axis

```r
pca <- prcomp(pc_in[cc, ], scale. = TRUE)
dp$PC1 <- pca$x[,1]; dp$PC2 <- pca$x[,2]; dp$PC3 <- pca$x[,3]
```

`prcomp(..., scale. = TRUE)` standardizes then decomposes the seven indices
(Eq. 4). `pca$x` holds the component **scores** (one value per patient per
component); columns 1–3 are attached to the data so PC1/PC2/PC3 can each be tested
as a predictor. Testing PC2 (the shape axis) and PC3 — which the manuscript never
did — is what shows the severity signal lives on the **size** axis (PC1), not the
shape axis.

### Block 7 — a validated severity index

```r
auc0 <- as.numeric(pROC::auc(dd$HH_hi, predict(g, type = "response"), quiet = TRUE))
opt <- replicate(500, {
  b  <- dd[sample(nrow(dd), replace = TRUE), ]
  gb <- glm(HH_hi ~ z_size_ratio + z_aspect_ratio + bleb, data = b, family = binomial)
  a_boot <- as.numeric(pROC::auc(b$HH_hi,  predict(gb, type="response"), quiet = TRUE))
  a_orig <- as.numeric(pROC::auc(dd$HH_hi, predict(gb, newdata = dd, type="response"), quiet = TRUE))
  a_boot - a_orig
})
cat(sprintf("Optimism-corrected AUC = %.3f\n", auc0 - mean(opt)))
```

`auc0` is the **apparent** discrimination — optimistic because the model was
scored on the same data it was fit on. The `replicate(500, …)` loop is the
Harrell bootstrap (Eq. 7): each iteration draws a bootstrap sample
`dd[sample(..., replace = TRUE), ]`, refits, and measures how much better the
model looks on its own bootstrap sample (`a_boot`) than on the original data
(`a_orig`). The average gap is the **optimism**; subtracting it gives an
honest, internally-validated AUC.
