## 02_analysis.R
## Corrected + extended analysis. Addresses reviewer critiques and explores
## stronger analyses. Predictors are standardized (per-SD) so ORs are comparable.
## Ordinal outcomes -> proportional-odds models (not linear).

suppressMessages({
  library(MASS); library(ordinal); library(mice); library(brant); library(pROC)
})
options(width = 120)
d <- readRDS("../outputs/analysis_data.rds")

COVARS <- c("age", "sex", "smoking_ever", "htn", "dm", "statin", "circulation")
INDICES <- c("z_size_ratio", "z_aspect_ratio", "z_maxdiam", "z_log_volume",
             "z_neck", "z_hw_ratio", "z_neck_parent")
OUTS <- list(mFS = "mFS_ord", HH = "HH_ord", WFNS = "WFNS_ord")

cat("\n########## 1. ADJUSTED PROPORTIONAL-ODDS ORDINAL REGRESSION (complete case) ##########\n")
cat("   (per-SD standardized predictor; OR>1 = higher severity; addresses Reviewer #2)\n\n")
res <- data.frame()
for (on in names(OUTS)) {
  for (pr in INDICES) {
    f <- as.formula(paste0(OUTS[[on]], " ~ ", pr, " + ", paste(COVARS, collapse = " + ")))
    m <- tryCatch(MASS::polr(f, data = d, Hess = TRUE, method = "logistic"),
                  error = function(e) NULL)
    if (is.null(m)) next
    ct <- coef(summary(m)); if (!pr %in% rownames(ct)) next
    b <- ct[pr, 1]; z <- ct[pr, 3]; p <- 2 * pnorm(abs(z), lower.tail = FALSE)
    res <- rbind(res, data.frame(outcome = on, predictor = sub("^z_", "", pr),
                 N = nrow(model.frame(m)), OR = exp(b),
                 lo = exp(b - 1.96 * ct[pr, 2]), hi = exp(b + 1.96 * ct[pr, 2]), p = p))
  }
}
res$p_FDR <- p.adjust(res$p, "BH")
res[, c("OR","lo","hi","p","p_FDR")] <- round(res[, c("OR","lo","hi","p","p_FDR")], 3)
print(res, row.names = FALSE)

cat("\n########## 2. SAME MODELS WITH MULTIPLE IMPUTATION (mice m=20) ##########\n")
cat("   (recovers the ~40% N lost to HTN missingness; addresses Reviewer #1 missingness)\n\n")
imp_vars <- unique(c(unlist(OUTS), INDICES, COVARS))
dm_imp <- d[, imp_vars]
# ordinal outcomes imputed as ordered factors
ini <- mice(dm_imp, maxit = 0, printFlag = FALSE)
meth <- ini$method
imp <- mice(dm_imp, m = 20, method = meth, printFlag = FALSE, seed = 2024)

## manual Rubin's rules pooling for one target coefficient across polr fits
pool_polr <- function(f, pr, imp) {
  ests <- c(); vars <- c(); dfres <- c()
  for (i in 1:imp$m) {
    di <- complete(imp, i)
    m <- tryCatch(MASS::polr(f, data = di, Hess = TRUE), error = function(e) NULL)
    if (is.null(m)) next
    ct <- coef(summary(m)); if (!pr %in% rownames(ct)) next
    ests <- c(ests, ct[pr, 1]); vars <- c(vars, ct[pr, 2]^2)
    dfres <- c(dfres, nrow(model.frame(m)))
  }
  if (length(ests) < 2) return(NULL)
  qbar <- mean(ests); ubar <- mean(vars); b <- var(ests)
  T <- ubar + (1 + 1/length(ests)) * b
  se <- sqrt(T); p <- 2 * pnorm(abs(qbar / se), lower.tail = FALSE)
  c(OR = exp(qbar), lo = exp(qbar - 1.96*se), hi = exp(qbar + 1.96*se), p = p)
}
res_mi <- data.frame()
for (on in names(OUTS)) {
  for (pr in c("z_size_ratio", "z_aspect_ratio", "z_maxdiam", "z_log_volume")) {
    f <- as.formula(paste0(OUTS[[on]], " ~ ", pr, " + ", paste(COVARS, collapse = " + ")))
    r <- pool_polr(f, pr, imp)
    if (is.null(r)) next
    res_mi <- rbind(res_mi, data.frame(outcome = on, predictor = sub("^z_", "", pr),
                 OR = r["OR"], lo = r["lo"], hi = r["hi"], p = r["p"]))
  }
}
res_mi$p_FDR <- p.adjust(res_mi$p, "BH")
res_mi[, c("OR","lo","hi","p","p_FDR")] <- round(res_mi[, c("OR","lo","hi","p","p_FDR")], 3)
print(res_mi, row.names = FALSE)

cat("\n########## 3. PROPORTIONAL-ODDS ASSUMPTION (Brant test, key models) ##########\n")
for (on in names(OUTS)) {
  f <- as.formula(paste0(OUTS[[on]], " ~ z_size_ratio + age + sex + htn"))
  m <- tryCatch(MASS::polr(f, data = d, Hess = TRUE), error = function(e) NULL)
  if (is.null(m)) next
  bt <- tryCatch(brant::brant(m), error = function(e) NULL)
  if (!is.null(bt)) cat(sprintf("  %s ~ size_ratio: Omnibus Brant p = %.3f (>0.05 = PO holds)\n",
                                on, bt["Omnibus", "probability"]))
}

cat("\n########## 4. PC1, PC2, PC3 vs OUTCOMES (paper only tested PC1) ##########\n")
pc_in <- d[, c("neck","maxdiam","volume","aspect_ratio","hw_ratio","size_ratio","neck_parent")]
cc <- complete.cases(pc_in)
pca <- prcomp(pc_in[cc, ], scale. = TRUE)
cat("Variance explained:", round(summary(pca)$importance[2, 1:4], 3), "\n")
cat("PC1 loadings:", paste(names(pc_in), round(pca$rotation[,1], 2), sep="=", collapse="  "), "\n")
cat("PC2 loadings:", paste(names(pc_in), round(pca$rotation[,2], 2), sep="=", collapse="  "), "\n")
dp <- d[cc, ]; dp$PC1 <- pca$x[,1]; dp$PC2 <- pca$x[,2]; dp$PC3 <- pca$x[,3]
for (on in names(OUTS)) {
  for (pc in c("PC1","PC2","PC3")) {
    f <- as.formula(paste0(OUTS[[on]], " ~ ", pc, " + age + sex + htn"))
    m <- tryCatch(MASS::polr(f, data = dp, Hess = TRUE), error = function(e) NULL)
    if (is.null(m)) next
    ct <- coef(summary(m)); p <- 2*pnorm(abs(ct[pc,3]), lower.tail=FALSE)
    cat(sprintf("  %-5s ~ %s: OR=%.2f p=%.3f\n", on, pc, exp(ct[pc,1]), p))
  }
}
# GCS (continuous) -- the one place composite fired in the paper
gc <- lm(GCS ~ PC1 + age + sex + htn, data = dp)
cat(sprintf("  GCS ~ PC1 (linear): beta=%.2f p=%.3f\n",
            coef(summary(gc))["PC1",1], coef(summary(gc))["PC1",4]))

cat("\n########## 5. BLEB / IRREGULARITY (in data, unused in paper) ##########\n")
d$bleb_f <- factor(ifelse(d$bleb == 1, "Yes", "No"))
for (on in names(OUTS)) {
  m <- tryCatch(MASS::polr(as.formula(paste0(OUTS[[on]], " ~ bleb_f + age + sex + htn")),
                data = d, Hess = TRUE), error = function(e) NULL)
  if (is.null(m)) next
  ct <- coef(summary(m)); rn <- grep("bleb", rownames(ct), value = TRUE)[1]
  p <- 2*pnorm(abs(ct[rn,3]), lower.tail=FALSE)
  cat(sprintf("  %-5s ~ bleb: OR=%.2f p=%.3f\n", on, exp(ct[rn,1]), p))
}
# size_ratio x bleb interaction (does irregularity amplify size effect?)
mi2 <- tryCatch(MASS::polr(HH_ord ~ z_size_ratio * bleb_f + age + sex, data = d, Hess = TRUE),
                error = function(e) NULL)
if (!is.null(mi2)) {
  ct <- coef(summary(mi2)); ix <- grep(":", rownames(ct), value = TRUE)[1]
  if (!is.na(ix)) cat(sprintf("  HH: size_ratio x bleb interaction p = %.3f\n",
                              2*pnorm(abs(ct[ix,3]), lower.tail=FALSE)))
}

cat("\n########## 6. NONLINEARITY / THRESHOLD for size ratio (HH) ##########\n")
d$SRq <- cut(d$size_ratio, breaks = quantile(d$size_ratio, c(0,.33,.67,1), na.rm=TRUE),
             labels = c("low","mid","high"), include.lowest = TRUE)
mt <- tryCatch(MASS::polr(HH_ord ~ SRq + age + sex + htn, data = d, Hess = TRUE),
               error = function(e) NULL)
if (!is.null(mt)) { ct <- coef(summary(mt))
  for (r in grep("SRq", rownames(ct), value = TRUE))
    cat(sprintf("  HH: %s vs low  OR=%.2f p=%.3f\n", r, exp(ct[r,1]),
                2*pnorm(abs(ct[r,3]), lower.tail=FALSE))) }

cat("\n########## 7. PARSIMONIOUS MULTIVARIABLE 'SEVERITY INDEX' + bootstrap c-index ##########\n")
cat("   (HH 3-5 vs 1-2 ~ size_ratio + aspect_ratio + bleb; internal-validated AUC)\n")
dd <- d[complete.cases(d[, c("HH_hi","z_size_ratio","z_aspect_ratio","bleb")]), ]
g <- glm(HH_hi ~ z_size_ratio + z_aspect_ratio + bleb, data = dd, family = binomial)
auc0 <- as.numeric(pROC::auc(dd$HH_hi, predict(g, type = "response"), quiet = TRUE))
set.seed(1)
opt <- replicate(500, {
  b <- dd[sample(nrow(dd), replace = TRUE), ]
  gb <- glm(HH_hi ~ z_size_ratio + z_aspect_ratio + bleb, data = b, family = binomial)
  a_boot <- as.numeric(pROC::auc(b$HH_hi, predict(gb, type="response"), quiet = TRUE))
  a_orig <- as.numeric(pROC::auc(dd$HH_hi, predict(gb, newdata = dd, type="response"), quiet = TRUE))
  a_boot - a_orig
})
cat(sprintf("  Apparent AUC=%.3f  Optimism=%.3f  Optimism-corrected AUC=%.3f  (N=%d)\n",
            auc0, mean(opt), auc0 - mean(opt), nrow(dd)))
print(round(summary(g)$coefficients, 3))

cat("\n########## DONE ##########\n")
