## 03_rwe_checks.R
## Real-world-evidence audit of the imputation + confounding assumptions.
## (1) Is MICE done correctly / does it inject bias?  (2) MNAR sensitivity.
## (3) E-values for unmeasured confounding.  (4) Missingness-mechanism diagnostic.

suppressMessages({ library(MASS); library(mice); library(EValue) })
options(width = 120)
d <- readRDS("../outputs/analysis_data.rds")

COVARS  <- c("age","sex","smoking_ever","htn","dm","statin","circulation")
OUTS    <- list(mFS="mFS_ord", HH="HH_ord", WFNS="WFNS_ord")
KEY     <- c("z_size_ratio","z_aspect_ratio")

## ---- helper: fit polr, return log-OR + SE for target term ---------------
fit1 <- function(f, pr, data) {
  m <- tryCatch(MASS::polr(f, data = data, Hess = TRUE), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  ct <- coef(summary(m)); if (!pr %in% rownames(ct)) return(NULL)
  c(b = ct[pr,1], se = ct[pr,2], n = nrow(model.frame(m)))
}

cat("############ AUDIT 1: Is the OUTCOME in the imputation model? ############\n")
imp_vars <- unique(c(unlist(OUTS), "GCS", "z_size_ratio","z_aspect_ratio",
                     "z_maxdiam","z_log_volume", COVARS))
cat("Imputation variable set includes outcomes:",
    all(unlist(OUTS) %in% imp_vars), "\n")
cat("--> Correct: including the outcome (and the correlated sister-scales +GCS as\n")
cat("    auxiliaries) is REQUIRED; omitting it biases associations toward the null.\n\n")

dm_imp <- d[, imp_vars]
ini  <- mice(dm_imp, maxit = 0, printFlag = FALSE)
imp  <- mice(dm_imp, m = 20, method = ini$method, printFlag = FALSE, seed = 2024)
cat("Imputation methods used per variable:\n"); print(ini$method[ini$method != ""])
cat("\nLogged events (collinearity/constant problems flagged by mice):\n")
print(if (is.null(imp$loggedEvents)) "none" else head(imp$loggedEvents, 10))

cat("\n############ AUDIT 2: Does MI inflate EFFECT SIZE, or only recover power? ############\n")
cat("If MI is unbiased, the pooled OR should be CLOSE to complete-case; only the CI\n")
cat("should tighten (more effective N). A large jump in the OR itself = red flag.\n\n")
cat(sprintf("%-6s %-14s %8s %8s %8s %6s %6s\n","out","pred","OR_CC","OR_MI","d_OR%","FMI","lambda"))
for (on in names(OUTS)) for (pr in KEY) {
  f <- as.formula(paste0(OUTS[[on]]," ~ ",pr," + ",paste(COVARS,collapse=" + ")))
  cc <- fit1(f, pr, d)
  ests <- c(); vars <- c()
  for (i in 1:imp$m) { r <- fit1(f, pr, complete(imp,i)); if(!is.null(r)){ ests<-c(ests,r["b"]); vars<-c(vars,r["se"]^2) } }
  m <- length(ests); qbar<-mean(ests); ubar<-mean(vars); B<-var(ests)
  Tt <- ubar + (1+1/m)*B; lambda <- (1+1/m)*B/Tt
  dfobs <- cc["n"]; fmi <- (lambda + 2/(dfobs+3))/(lambda+1)
  cat(sprintf("%-6s %-14s %8.2f %8.2f %+7.1f %6.2f %6.2f\n",
      on, sub("^z_","",pr), exp(cc["b"]), exp(qbar),
      100*(exp(qbar)-exp(cc["b"]))/exp(cc["b"]), fmi, lambda))
}

cat("\n############ AUDIT 3: MI-then-DELETE (von Hippel) sensitivity ############\n")
cat("Impute covariates using the outcome, but analyze only OBSERVED-outcome rows.\n")
cat("If conclusions match full-MI, the imputed-outcome rows are not driving them.\n\n")
for (on in names(OUTS)) for (pr in KEY) {
  f <- as.formula(paste0(OUTS[[on]]," ~ ",pr," + ",paste(COVARS,collapse=" + ")))
  ests<-c(); vars<-c(); ocol <- OUTS[[on]]
  for (i in 1:imp$m) {
    di <- complete(imp,i); di <- di[!is.na(d[[ocol]]), ]   # keep only truly-observed outcome
    r <- fit1(f, pr, di); if(!is.null(r)){ ests<-c(ests,r["b"]); vars<-c(vars,r["se"]^2) }
  }
  m<-length(ests); qbar<-mean(ests); Tt<-mean(vars)+(1+1/m)*var(ests); se<-sqrt(Tt)
  p<-2*pnorm(abs(qbar/se),lower.tail=FALSE)
  cat(sprintf("  MID  %-5s ~ %-12s OR=%.2f (%.2f-%.2f) p=%.3f\n",
      on, sub("^z_","",pr), exp(qbar), exp(qbar-1.96*se), exp(qbar+1.96*se), p))
}

cat("\n############ AUDIT 4: MNAR delta-adjustment (pattern-mixture) ############\n")
cat("Worst case: patients with MISSING grades are systematically SICKER. Shift every\n")
cat("imputed grade up by delta categories; re-test size_ratio -> HH. Does it survive?\n\n")
for (delta in c(0, 0.5, 1)) {
  ests<-c(); vars<-c()
  for (i in 1:imp$m) {
    di <- complete(imp,i)
    miss <- is.na(d$HH_ord)
    hh <- as.integer(di$HH_ord); hh[miss] <- pmin(5, hh[miss] + delta)
    di$HHc <- hh
    r <- fit1(HHc ~ z_size_ratio + age + sex + htn,
              "z_size_ratio", within(di, HHc <- factor(round(HHc), ordered=TRUE)))
    if(!is.null(r)){ ests<-c(ests,r["b"]); vars<-c(vars,r["se"]^2) }
  }
  m<-length(ests); qbar<-mean(ests); se<-sqrt(mean(vars)+(1+1/m)*var(ests))
  cat(sprintf("  delta=+%.1f grade: size_ratio->HH  OR=%.2f  p=%.3f\n",
      delta, exp(qbar), 2*pnorm(abs(qbar/se),lower.tail=FALSE)))
}

cat("\n############ RWE ADD-ON: E-VALUES for unmeasured confounding ############\n")
cat("How strong (on the OR scale) must an unmeasured confounder be, associated with\n")
cat("BOTH size ratio and severity, to explain away the association? (outcomes common\n")
cat("-> rare=FALSE). Directly answers the 'unmeasured confounding' limitation.\n\n")
mi_or <- function(on, pr) {
  f <- as.formula(paste0(OUTS[[on]]," ~ ",pr," + ",paste(COVARS,collapse=" + ")))
  ests<-c(); vars<-c()
  for (i in 1:imp$m){ r<-fit1(f,pr,complete(imp,i)); if(!is.null(r)){ests<-c(ests,r["b"]);vars<-c(vars,r["se"]^2)} }
  m<-length(ests); qbar<-mean(ests); se<-sqrt(mean(vars)+(1+1/m)*var(ests))
  c(or=exp(qbar), lo=exp(qbar-1.96*se), hi=exp(qbar+1.96*se))
}
for (on in c("mFS","HH","WFNS")) for (pr in KEY) {
  e <- mi_or(on, pr)
  ev <- EValue::evalues.OR(est = e["or"], lo = e["lo"], hi = e["hi"], rare = FALSE)
  cat(sprintf("  %-5s ~ %-12s OR=%.2f (%.2f-%.2f)  E-value point=%.2f  CI=%.2f\n",
      on, sub("^z_","",pr), e["or"], e["lo"], e["hi"],
      ev["E-values","point"], ev["E-values","lower"]))
}

cat("\n############ AUDIT 5: Missingness mechanism (MCAR vs MAR diagnostic) ############\n")
cat("Does observed morphology predict WHETHER the grade is missing? If yes, data are\n")
cat("MAR-recoverable (we condition on morphology); if nothing predicts it, ~MCAR.\n\n")
d$HH_missing <- as.integer(is.na(d$HH))
gm <- glm(HH_missing ~ z_size_ratio + z_maxdiam + age + sex, data = d, family = binomial)
print(round(summary(gm)$coefficients, 3))

cat("\n############ DONE ############\n")
