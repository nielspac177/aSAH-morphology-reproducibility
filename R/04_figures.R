## 04_figures.R
## Publication figures (300 dpi TIFF + vector PDF), colorblind-safe (Okabe-Ito).
## Fig 1  Forest plot of adjusted per-SD ordinal ORs, complete-case vs MICE.
## Fig 2  PCA: size axis (PC1) vs shape axis (PC2) + scree + PC->outcome.
## Fig 3  Sensitivity: E-values, size-ratio threshold, MNAR delta.

suppressMessages({
  library(MASS); library(mice); library(EValue)
  library(ggplot2); library(patchwork); library(ggrepel); library(dplyr); library(tidyr)
})
d <- readRDS("../outputs/analysis_data.rds")
dir.create("../figures", showWarnings = FALSE)

OK <- c(blue="#0072B2", orange="#E69F00", green="#009E73", red="#D55E00",
        purple="#CC79A7", sky="#56B4E9", grey="#999999", black="#000000")
theme_pub <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom")
theme_set(theme_pub)

COVARS <- c("age","sex","smoking_ever","htn","dm","statin","circulation")
INDICES <- c(size_ratio="z_size_ratio", aspect_ratio="z_aspect_ratio",
             max_diameter="z_maxdiam", volume="z_log_volume",
             neck_width="z_neck", hw_ratio="z_hw_ratio", neck_parent="z_neck_parent")
OUTS <- c(mFS="mFS_ord", HH="HH_ord", WFNS="WFNS_ord")

fit1 <- function(f, pr, data) {
  m <- tryCatch(MASS::polr(f, data = data, Hess = TRUE), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  ct <- coef(summary(m)); if (!pr %in% rownames(ct)) return(NULL)
  c(b = ct[pr,1], se = ct[pr,2])
}

## ---------- compute complete-case + MICE ORs -----------------------------
message("Fitting complete-case models ...")
cc <- expand.grid(index = names(INDICES), outcome = names(OUTS),
                  stringsAsFactors = FALSE)
cc$OR <- NA; cc$lo <- NA; cc$hi <- NA; cc$p <- NA
for (i in seq_len(nrow(cc))) {
  pr <- INDICES[[cc$index[i]]]
  f <- as.formula(paste0(OUTS[[cc$outcome[i]]], " ~ ", pr, " + ", paste(COVARS, collapse=" + ")))
  r <- fit1(f, pr, d); if (is.null(r)) next
  cc$OR[i] <- exp(r["b"]); cc$lo[i] <- exp(r["b"]-1.96*r["se"])
  cc$hi[i] <- exp(r["b"]+1.96*r["se"]); cc$p[i] <- 2*pnorm(abs(r["b"]/r["se"]), lower.tail=FALSE)
}
cc$method <- "Complete-case"

message("Multiple imputation (m=20) ...")
imp_vars <- unique(c(unname(OUTS), "GCS", unname(INDICES), COVARS))
ini <- mice(d[, imp_vars], maxit = 0, printFlag = FALSE)
imp <- mice(d[, imp_vars], m = 20, method = ini$method, printFlag = FALSE, seed = 2024)
mi <- cc[, c("index","outcome")]; mi$OR<-NA; mi$lo<-NA; mi$hi<-NA; mi$p<-NA
for (i in seq_len(nrow(mi))) {
  pr <- INDICES[[mi$index[i]]]
  f <- as.formula(paste0(OUTS[[mi$outcome[i]]], " ~ ", pr, " + ", paste(COVARS, collapse=" + ")))
  ests<-c(); vars<-c()
  for (k in 1:imp$m){ r<-fit1(f,pr,complete(imp,k)); if(!is.null(r)){ests<-c(ests,r["b"]);vars<-c(vars,r["se"]^2)} }
  m<-length(ests); qbar<-mean(ests); se<-sqrt(mean(vars)+(1+1/m)*var(ests))
  mi$OR[i]<-exp(qbar); mi$lo[i]<-exp(qbar-1.96*se); mi$hi[i]<-exp(qbar+1.96*se)
  mi$p[i]<-2*pnorm(abs(qbar/se), lower.tail=FALSE)
}
mi$method <- "Multiple imputation"

fp <- bind_rows(cc, mi)
fp$p_FDR <- ave(fp$p, fp$method, FUN = function(p) p.adjust(p, "BH"))
fp$sig <- ifelse(fp$p_FDR < 0.05, "FDR < 0.05", ifelse(fp$p < 0.05, "p < 0.05", "NS"))
lab <- c(size_ratio="Size ratio", aspect_ratio="Aspect ratio", max_diameter="Max diameter",
         volume="Volume (log)", neck_width="Neck width", hw_ratio="Height/width ratio",
         neck_parent="Neck/parent ratio")
fp$index <- factor(lab[fp$index], levels = rev(lab))
fp$outcome <- factor(fp$outcome, levels = c("mFS","HH","WFNS"),
                     labels = c("modified Fisher","Hunt & Hess","WFNS"))
fp$method <- factor(fp$method, levels = c("Complete-case","Multiple imputation"))

## ============ FIGURE 1 : JAMA/NEJM-style table with embedded forest =======
## Primary (multiple-imputation) estimates. Left = text columns; right = forest,
## row-aligned; grouped by outcome. Complete-case ORs shown as a text column.
ind_order <- c("Size ratio","Aspect ratio","Max diameter","Volume (log)",
               "Neck/parent ratio","Neck width","Height/width ratio")
out_order <- c("modified Fisher","Hunt & Hess","WFNS")
mi_ <- subset(fp, method == "Multiple imputation")
cc_ <- subset(fp, method == "Complete-case",
              select = c(index, outcome, OR, lo, hi))
names(cc_)[3:5] <- c("ccOR","cclo","cchi")
tab <- merge(mi_, cc_, by = c("index","outcome"))
tab$index   <- factor(as.character(tab$index),   levels = ind_order)
tab$outcome <- factor(as.character(tab$outcome), levels = out_order)
tab <- tab[order(tab$outcome, tab$index), ]

# assign y positions: a header row per outcome group, then its 7 index rows
rows <- data.frame(); y <- 0
for (o in out_order) {
  y <- y - 1
  rows <- rbind(rows, data.frame(y = y, kind = "header", outcome = o,
    index = NA, txt_pred = as.character(o), OR = NA, lo = NA, hi = NA,
    sig = NA, ci = NA, pval = NA, fdrv = NA))
  sub <- tab[tab$outcome == o, ]
  for (i in ind_order) {
    r <- sub[sub$index == i, ]; y <- y - 1
    rows <- rbind(rows, data.frame(y = y, kind = "data", outcome = o,
      index = i, txt_pred = paste0("   ", i), OR = r$OR, lo = r$lo, hi = r$hi,
      sig = r$sig,
      ci = sprintf("%.2f (%.2f–%.2f)", r$OR, r$lo, r$hi),
      pval = ifelse(r$p < .001, "<0.001", sprintf("%.3f", r$p)),
      fdrv = ifelse(r$p_FDR < .001, "<0.001", sprintf("%.3f", r$p_FDR))))
  }
}
ytop <- 1
colcol <- c("FDR < 0.05"=unname(OK["red"]), "p < 0.05"=unname(OK["blue"]), "NS"=unname(OK["grey"]))

# ---- left panel: the table ----
dat <- subset(rows, kind == "data"); hdr <- subset(rows, kind == "header")
tbl <- ggplot() +
  # column headers
  annotate("text", x = 0,    y = ytop, label = "Predictor",     hjust = 0, fontface = "bold", size = 3) +
  annotate("text", x = 3.05, y = ytop, label = "OR (95% CI)",   hjust = 0, fontface = "bold", size = 3) +
  annotate("text", x = 5.5,  y = ytop, label = "p",             hjust = 0, fontface = "bold", size = 3) +
  annotate("text", x = 6.35, y = ytop, label = "FDR",           hjust = 0, fontface = "bold", size = 3) +
  geom_text(data = hdr, aes(0, y, label = txt_pred), hjust = 0, fontface = "bold", size = 3) +
  geom_text(data = dat, aes(0, y, label = txt_pred), hjust = 0, size = 2.9) +
  geom_text(data = dat, aes(3.05, y, label = ci),   hjust = 0, size = 2.9) +
  geom_text(data = dat, aes(5.5,  y, label = pval), hjust = 0, size = 2.9) +
  geom_text(data = dat, aes(6.35, y, label = fdrv, color = sig), hjust = 0, size = 2.9, show.legend = FALSE) +
  scale_color_manual(values = colcol) +
  coord_cartesian(xlim = c(0, 7.4), ylim = c(min(rows$y) - 0.6, ytop + 0.6), expand = FALSE) +
  theme_void()

# ---- right panel: the forest ----
brks <- c(0.5, 1, 2, 4)
fst <- ggplot(dat, aes(OR, y, color = sig)) +
  annotate("text", x = 1, y = ytop, label = "Odds ratio (per 1 SD)", fontface = "bold", size = 3) +
  geom_vline(xintercept = 1, linetype = 2, color = OK["grey"]) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.6) +
  geom_point(size = 2.2) +
  geom_hline(data = hdr, aes(yintercept = y + 0.5), color = "grey85", linewidth = 0.3) +
  scale_x_log10(breaks = brks, labels = brks) +
  scale_color_manual(values = colcol, name = NULL) +
  coord_cartesian(xlim = c(0.5, 4), ylim = c(min(rows$y) - 0.6, ytop + 0.6), expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        axis.text.y = element_blank(), legend.position = "bottom")

f1 <- tbl + fst + plot_layout(widths = c(1.15, 1)) +
  plot_annotation(
    title = "Adjusted association of aneurysm morphology with aSAH severity",
    subtitle = "Proportional-odds ordinal regression, multiple imputation; OR > 1 = higher severity. Complete-case estimates in Table 2.",
    theme = theme(plot.title = element_text(face = "bold", size = 12)))
ggsave("../figures/Figure1_forest.tiff", f1, width = 9.5, height = 7.2, dpi = 300, compression = "lzw", type = "cairo")
ggsave("../figures/Figure1_forest.pdf",  f1, width = 9.5, height = 7.2)
ggsave("../figures/Figure1_forest.png",  f1, width = 9.5, height = 7.2, dpi = 150, type = "cairo")  # for README

## ================= FIGURE 2 : PCA size vs shape ==========================
pc_in <- d[, c("neck","maxdiam","volume","aspect_ratio","hw_ratio","size_ratio","neck_parent")]
cc_pc <- complete.cases(pc_in)
pca <- prcomp(pc_in[cc_pc, ], scale. = TRUE)
ve <- summary(pca)$importance[2, 1:7]
load <- as.data.frame(pca$rotation[, 1:2]); load$var <- c("Neck width","Max diameter","Volume",
  "Aspect ratio","Height/width","Size ratio","Neck/parent")
p2a <- ggplot(load, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, color = OK["grey"], linetype = 2) +
  geom_vline(xintercept = 0, color = OK["grey"], linetype = 2) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.18,"cm")), color = OK["blue"]) +
  geom_text_repel(aes(label = var), size = 3.2, seg.color = "grey70") +
  labs(title = "A  Morphology axes",
       x = sprintf("PC1 — size axis (%.0f%%)", 100*ve[1]),
       y = sprintf("PC2 — shape axis (%.0f%%)", 100*ve[2]))
scree <- data.frame(PC = factor(paste0("PC",1:7), levels=paste0("PC",1:7)), ve = as.numeric(ve))
p2b <- ggplot(scree, aes(PC, ve)) +
  geom_col(fill = OK["sky"]) +
  geom_text(aes(label = sprintf("%.0f%%",100*ve)), vjust = -0.4, size = 3) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult=c(0,0.12))) +
  labs(title = "B  Variance explained", x = NULL, y = "Proportion of variance")
# PC1/PC2 -> outcomes
dp <- d[cc_pc, ]; dp$PC1 <- pca$x[,1]; dp$PC2 <- pca$x[,2]
pcres <- data.frame()
for (on in names(OUTS)) for (pc in c("PC1","PC2")) {
  f <- as.formula(paste0(OUTS[[on]], " ~ ", pc, " + age + sex + htn"))
  r <- fit1(f, pc, dp); if (is.null(r)) next
  pcres <- rbind(pcres, data.frame(outcome=on, comp=pc, OR=exp(r["b"]),
    lo=exp(r["b"]-1.96*r["se"]), hi=exp(r["b"]+1.96*r["se"]),
    p=2*pnorm(abs(r["b"]/r["se"]), lower.tail=FALSE)))
}
pcres$outcome <- factor(pcres$outcome, levels=c("mFS","HH","WFNS"))
pcres$comp <- factor(pcres$comp, labels=c("PC1 (size)","PC2 (shape)"))
p2c <- ggplot(pcres, aes(OR, outcome, color = comp)) +
  geom_vline(xintercept = 1, linetype = 2, color = OK["grey"]) +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.2, position=position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5), size = 2) +
  scale_color_manual(values = c(unname(OK["red"]), unname(OK["blue"])), name=NULL) +
  scale_x_log10() +
  labs(title = "C  Component → severity", x = "OR per 1 SD (log scale)", y = NULL)
f2 <- (p2a | p2b) / p2c + plot_layout(heights = c(1.3, 1))
ggsave("../figures/Figure2_PCA.tiff", f2, width = 9, height = 8, dpi = 300, compression = "lzw", type = "cairo")
ggsave("../figures/Figure2_PCA.pdf",  f2, width = 9, height = 8)

## ================= FIGURE 3 : sensitivity ================================
# (a) E-values for the two primary indices (from MI ORs)
prim <- fp %>% filter(method=="Multiple imputation",
                      index %in% c("Size ratio","Aspect ratio"))
ev <- do.call(rbind, lapply(seq_len(nrow(prim)), function(i){
  e <- EValue::evalues.OR(prim$OR[i], prim$lo[i], prim$hi[i], rare = FALSE)
  data.frame(label = paste(prim$outcome[i], prim$index[i], sep=" / "),
             point = e["E-values","point"], ci = e["E-values","lower"])
}))
ev$ci[is.na(ev$ci)] <- 1
evl <- pivot_longer(ev, c(point, ci), names_to = "type", values_to = "E")
evl$type <- factor(evl$type, levels=c("point","ci"), labels=c("Point estimate","CI limit"))
p3a <- ggplot(evl, aes(E, reorder(label, E), fill = type)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_vline(xintercept = 1, linetype = 2, color = OK["grey"]) +
  scale_fill_manual(values = c(unname(OK["orange"]), unname(OK["sky"])), name=NULL) +
  labs(title = "A  E-values (unmeasured confounding)",
       x = "Confounder association needed (OR)", y = NULL)
# (b) size-ratio tertiles -> HH probability of high grade
d$SRq <- cut(d$size_ratio, breaks = quantile(d$size_ratio, c(0,.33,.67,1), na.rm=TRUE),
             labels = c("Low","Mid","High"), include.lowest = TRUE)
thr <- d %>% filter(!is.na(SRq), !is.na(HH)) %>% group_by(SRq) %>%
  summarise(p_hi = mean(HH>=3), n = n(),
            se = sqrt(p_hi*(1-p_hi)/n), .groups="drop")
p3b <- ggplot(thr, aes(SRq, p_hi)) +
  geom_col(fill = OK["green"], width = 0.65) +
  geom_errorbar(aes(ymin = pmax(0,p_hi-1.96*se), ymax = pmin(1,p_hi+1.96*se)), width = 0.2) +
  geom_text(aes(label = sprintf("%.0f%%\n(n=%d)",100*p_hi,n)), vjust = -0.3, size = 3) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  labs(title = "B  Size-ratio threshold", x = "Size-ratio tertile",
       y = "High-grade HH (3–5)")
# (c) MNAR delta robustness for size-ratio -> HH
deltas <- c(0, 0.5, 1, 1.5)
mn <- data.frame()
for (dl in deltas) {
  ests<-c(); vars<-c()
  for (k in 1:imp$m) {
    di <- complete(imp,k); miss <- is.na(d$HH_ord)
    hh <- as.integer(di$HH_ord); hh[miss] <- pmin(5, hh[miss]+dl)
    di$HHc <- factor(round(hh), ordered=TRUE)
    r <- fit1(HHc ~ z_size_ratio + age + sex + htn, "z_size_ratio", di)
    if(!is.null(r)){ ests<-c(ests,r["b"]); vars<-c(vars,r["se"]^2) }
  }
  m<-length(ests); qbar<-mean(ests); se<-sqrt(mean(vars)+(1+1/m)*var(ests))
  mn <- rbind(mn, data.frame(delta=dl, OR=exp(qbar), lo=exp(qbar-1.96*se), hi=exp(qbar+1.96*se)))
}
p3c <- ggplot(mn, aes(factor(delta), OR)) +
  geom_hline(yintercept = 1, linetype = 2, color = OK["grey"]) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.15, color = OK["purple"]) +
  geom_point(size = 2.4, color = OK["purple"]) +
  labs(title = "C  MNAR sensitivity (size ratio → HH)",
       x = "Assumed severity shift for missing grades (categories)",
       y = "OR per 1 SD")
f3 <- p3a / (p3b | p3c) + plot_layout(heights = c(1, 1))
ggsave("../figures/Figure3_sensitivity.tiff", f3, width = 9, height = 8, dpi = 300, compression = "lzw", type = "cairo")
ggsave("../figures/Figure3_sensitivity.pdf",  f3, width = 9, height = 8)

message("Wrote Figure1_forest, Figure2_PCA, Figure3_sensitivity (.tiff + .pdf) to ../figures/")
