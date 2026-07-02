## 05_tables.R
## Publication tables, all regenerated from data. Writes CSVs to ../outputs/tables/
## and a combined Markdown file that is converted to Word with pandoc.

suppressMessages({ library(MASS); library(mice); library(EValue) })
options(width = 140)
d <- readRDS("../outputs/analysis_data.rds")
td <- "../outputs/tables"; dir.create(td, showWarnings = FALSE, recursive = TRUE)
MD <- c()                                   # markdown accumulator
add <- function(...) MD <<- c(MD, sprintf(...))

COVARS  <- c("age","sex","smoking_ever","htn","dm","statin","circulation")
INDICES <- c(`Size ratio`="z_size_ratio", `Aspect ratio`="z_aspect_ratio",
             `Max diameter`="z_maxdiam", `Volume (log)`="z_log_volume",
             `Neck width`="z_neck", `Height/width ratio`="z_hw_ratio",
             `Neck/parent ratio`="z_neck_parent")
OUTS <- c(`modified Fisher`="mFS_ord", `Hunt & Hess`="HH_ord", WFNS="WFNS_ord")

fit1 <- function(f, pr, data) {
  m <- tryCatch(MASS::polr(f, data = data, Hess = TRUE), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  ct <- coef(summary(m)); if (!pr %in% rownames(ct)) return(NULL)
  c(b = ct[pr,1], se = ct[pr,2], n = nrow(model.frame(m)))
}
fmt <- function(or, lo, hi, p) sprintf("%.2f (%.2f–%.2f), %s", or, lo, hi,
                                       ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

## ---- MICE ---------------------------------------------------------------
message("Imputing (m=20)...")
imp_vars <- unique(c(unname(OUTS), "GCS", unname(INDICES), COVARS))
ini <- mice(d[, imp_vars], maxit = 0, printFlag = FALSE)
imp <- mice(d[, imp_vars], m = 20, method = ini$method, printFlag = FALSE, seed = 2024)
mi_est <- function(f, pr) {
  ests <- c(); vars <- c()
  for (k in 1:imp$m){ r <- fit1(f, pr, complete(imp,k)); if(!is.null(r)){ests<-c(ests,r["b"]);vars<-c(vars,r["se"]^2)} }
  m<-length(ests); qbar<-mean(ests); se<-sqrt(mean(vars)+(1+1/m)*var(ests))
  c(or=exp(qbar), lo=exp(qbar-1.96*se), hi=exp(qbar+1.96*se), p=2*pnorm(abs(qbar/se), lower.tail=FALSE))
}

## ======================= TABLE 1 — descriptives ==========================
t1 <- data.frame(Variable=character(), Value=character())
push <- function(v, s) t1[nrow(t1)+1, ] <<- c(v, s)
mm <- function(x) sprintf("%.1f (%.1f)", mean(x,na.rm=TRUE), sd(x,na.rm=TRUE))
np <- function(x, lvl) { n<-sum(x==lvl,na.rm=TRUE); sprintf("%d (%.0f%%)", n, 100*n/sum(!is.na(x))) }
push("N", as.character(nrow(d)))
push("Age, years — mean (SD)", mm(d$age))
push("Female sex", np(d$sex,"Female"))
push("Current/former smoker", sprintf("%d (%.0f%%)", sum(d$smoking_ever=="Ever",na.rm=TRUE),
     100*mean(d$smoking_ever=="Ever",na.rm=TRUE)))
push("Hypertension", np(d$htn,"Yes")); push("Diabetes", np(d$dm,"Yes"))
push("Statin use", np(d$statin,"Yes")); push("Posterior circulation", np(d$circulation,"Posterior"))
push("HH 3–5", sprintf("%d (%.0f%%)", sum(d$HH>=3,na.rm=TRUE), 100*mean(d$HH>=3,na.rm=TRUE)))
push("mFS 3–4", sprintf("%d (%.0f%%)", sum(d$mFS>=3,na.rm=TRUE), 100*mean(d$mFS>=3,na.rm=TRUE)))
push("WFNS 4–5", sprintf("%d (%.0f%%)", sum(d$WFNS>=4,na.rm=TRUE), 100*mean(d$WFNS>=4,na.rm=TRUE)))
push("Max diameter, mm — mean (SD)", mm(d$maxdiam))
push("Aspect ratio — mean (SD)", sprintf("%.2f (%.2f)", mean(d$aspect_ratio,na.rm=TRUE), sd(d$aspect_ratio,na.rm=TRUE)))
push("Size ratio — mean (SD)", sprintf("%.2f (%.2f)", mean(d$size_ratio,na.rm=TRUE), sd(d$size_ratio,na.rm=TRUE)))
push("Volume, mm³ — mean (SD)", sprintf("%.0f (%.0f)", mean(d$volume,na.rm=TRUE), sd(d$volume,na.rm=TRUE)))
push("Bleb present", np(factor(ifelse(d$bleb==1,"Yes","No")),"Yes"))
write.csv(t1, file.path(td,"Table1_descriptives.csv"), row.names=FALSE)
add("## Table 1. Cohort characteristics (N = %d)\n", nrow(d))
add("| %s | %s |", names(t1)[1], names(t1)[2]); add("|---|---|")
for(i in 1:nrow(t1)) add("| %s | %s |", t1$Variable[i], t1$Value[i]); add("")

## ============ TABLE 2 — adjusted ordinal ORs (CC + MICE + FDR) ============
grid <- expand.grid(index=names(INDICES), outcome=names(OUTS), stringsAsFactors=FALSE)
res <- data.frame()
for (i in seq_len(nrow(grid))) {
  pr <- INDICES[[grid$index[i]]]
  f <- as.formula(paste0(OUTS[[grid$outcome[i]]], " ~ ", pr, " + ", paste(COVARS, collapse=" + ")))
  cc <- fit1(f, pr, d); mi <- mi_est(f, pr)
  cc_or <- exp(cc["b"]); cc_lo <- exp(cc["b"]-1.96*cc["se"]); cc_hi <- exp(cc["b"]+1.96*cc["se"])
  cc_p <- 2*pnorm(abs(cc["b"]/cc["se"]), lower.tail=FALSE)
  res <- rbind(res, data.frame(outcome=grid$outcome[i], index=grid$index[i], N=cc["n"],
    cc=fmt(cc_or,cc_lo,cc_hi,cc_p), cc_p=cc_p,
    mi=fmt(mi["or"],mi["lo"],mi["hi"],mi["p"]), mi_p=mi["p"]))
}
res$mi_fdr <- p.adjust(res$mi_p, "BH")
write.csv(res, file.path(td,"Table2_ordinal_ORs.csv"), row.names=FALSE)
add("## Table 2. Adjusted per-SD ordinal odds ratios for morphology → severity")
add("Proportional-odds ordinal regression, one index per model, adjusted for age, sex, smoking, hypertension, diabetes, statin, circulation. OR > 1 = higher severity.\n")
add("| Outcome | Predictor | N | Complete-case OR (95%% CI), p | Multiple-imputation OR (95%% CI), p | MI FDR |")
add("|---|---|---|---|---|---|")
for(i in 1:nrow(res)) add("| %s | %s | %d | %s | %s | %s |",
  res$outcome[i], res$index[i], res$N[i], res$cc[i], res$mi[i],
  ifelse(res$mi_fdr[i]<0.001,"<0.001",sprintf("%.3f",res$mi_fdr[i])))
add("")

## ============ TABLE 3 — FDR-significant under MICE (headline) =============
sig <- res[res$mi_fdr < 0.05, ]
sig <- sig[order(sig$mi_fdr), ]
write.csv(sig, file.path(td,"Table3_FDR_significant.csv"), row.names=FALSE)
add("## Table 3. Morphometric associations surviving FDR correction (multiple imputation)")
add("| Outcome | Predictor | MI OR (95%% CI), p | MI FDR |")
add("|---|---|---|---|")
for(i in 1:nrow(sig)) add("| %s | %s | %s | %s |", sig$outcome[i], sig$index[i], sig$mi[i],
  ifelse(sig$mi_fdr[i]<0.001,"<0.001",sprintf("%.3f",sig$mi_fdr[i])))
add("")

## ============ SUPP S1 — PC1/PC2/PC3 → outcomes (ordinal) ==================
pc_in <- d[, c("neck","maxdiam","volume","aspect_ratio","hw_ratio","size_ratio","neck_parent")]
ccpc <- complete.cases(pc_in); pca <- prcomp(pc_in[ccpc,], scale.=TRUE)
dp <- d[ccpc,]; dp$PC1<-pca$x[,1]; dp$PC2<-pca$x[,2]; dp$PC3<-pca$x[,3]
s1 <- data.frame()
for (on in names(OUTS)) for (pc in c("PC1","PC2","PC3")) {
  f <- as.formula(paste0(OUTS[[on]], " ~ ", pc, " + age + sex + htn"))
  r <- fit1(f, pc, dp); if(is.null(r)) next
  s1 <- rbind(s1, data.frame(outcome=on, component=pc,
    est=fmt(exp(r["b"]),exp(r["b"]-1.96*r["se"]),exp(r["b"]+1.96*r["se"]),2*pnorm(abs(r["b"]/r["se"]),lower.tail=FALSE))))
}
gc <- lm(GCS ~ PC1 + age + sex + htn, data=dp); gcc <- coef(summary(gc))["PC1",]
write.csv(s1, file.path(td,"SuppS1_PCA_components.csv"), row.names=FALSE)
add("## Supplementary Table S1. Principal components → severity (ordinal, adjusted)")
add("PC1 = size axis (48%% var); PC2 = shape axis (22%%).\n")
add("| Outcome | Component | OR (95%% CI), p |"); add("|---|---|---|")
for(i in 1:nrow(s1)) add("| %s | %s | %s |", s1$outcome[i], s1$component[i], s1$est[i])
add("| GCS (linear) | PC1 | β=%.2f, p=%.3f |", gcc[1], gcc[4]); add("")

## ============ SUPP S2 — E-values (primary indices, MICE) =================
add("## Supplementary Table S2. E-values for unmeasured confounding (MI estimates)")
add("| Outcome | Predictor | OR (95%% CI) | E-value (point) | E-value (CI) |"); add("|---|---|---|---|---|")
for (on in names(OUTS)) for (idx in c("Size ratio","Aspect ratio")) {
  pr <- INDICES[[idx]]
  f <- as.formula(paste0(OUTS[[on]], " ~ ", pr, " + ", paste(COVARS, collapse=" + ")))
  mi <- mi_est(f, pr); ev <- EValue::evalues.OR(mi["or"], mi["lo"], mi["hi"], rare=FALSE)
  ci <- ev["E-values","lower"]; ci <- ifelse(is.na(ci),1,ci)
  add("| %s | %s | %.2f (%.2f–%.2f) | %.2f | %.2f |", on, idx, mi["or"],mi["lo"],mi["hi"],
      ev["E-values","point"], ci)
}
add("")

writeLines(MD, file.path(td, "TABLES.md"))
message("Wrote CSVs + TABLES.md to ", td)

## ---- convert to Word ----------------------------------------------------
ok <- suppressWarnings(system2("pandoc",
  c(shQuote(file.path(td,"TABLES.md")), "-o", shQuote(file.path(td,"TABLES.docx"))),
  stdout=TRUE, stderr=TRUE))
if (file.exists(file.path(td,"TABLES.docx"))) message("Wrote TABLES.docx") else message("pandoc note: ", paste(ok, collapse=" "))
