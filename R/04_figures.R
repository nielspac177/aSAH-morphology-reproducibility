## 04_figures.R
## Publication figures — JAMA house style (palette/typeface from _journal_theme.R).
## Fig 1  Table-with-forest of adjusted per-SD ordinal ORs (multiple imputation).
## Fig 2  PCA: size axis (PC1) vs shape axis (PC2) + scree + PC->outcome.
## Fig 3  Sensitivity: E-values (dumbbell), size-ratio threshold, MNAR delta.

suppressMessages({
  library(MASS); library(mice); library(EValue)
  library(ggplot2); library(patchwork); library(ggrepel); library(dplyr); library(tidyr)
})
source("_journal_theme.R")
d <- readRDS("../outputs/analysis_data.rds")
dir.create("../figures", showWarnings = FALSE)

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

## ============ FIGURE 1 : JAMA table-with-forest ==========================
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

## ---- classic JAMA forest: one ggplot, monochrome, forest in the middle,
##      OR (95% CI) + FDR columns at right, grey row-group bands, favors axis.
INKF <- "#1A1A1A"                      # forest ink (near-black)
## map OR (log scale) to plot x within the forest column [fx0, fx1]
fx0 <- 30; fx1 <- 64; orlo <- 0.5; orhi <- 4
px <- function(or) fx0 + (log10(or) - log10(orlo)) / (log10(orhi) - log10(orlo)) * (fx1 - fx0)

recs <- list(); y <- 0; bands <- list()
for (o in out_order) {
  y <- y - 1.15
  bands[[length(bands)+1]] <- data.frame(y = y, outcome = o)
  sub <- tab[tab$outcome == o, ]
  for (i in ind_order) {
    r <- sub[sub$index == i, ]; y <- y - 1
    recs[[length(recs)+1]] <- data.frame(y = y, index = i,
      OR = r$OR, lo = r$lo, hi = r$hi,
      ci = sprintf("%.2f (%.2f–%.2f)", r$OR, r$lo, r$hi),
      fdrv = ifelse(r$p_FDR < .001, "<0.001", sprintf("%.3f", r$p_FDR)),
      fdr_sig = r$p_FDR < 0.05,
      wt = 1 / (log(r$hi) - log(r$lo)))            # precision -> marker size
  }
  y <- y - 0.55
}
dat <- do.call(rbind, recs); bnd <- do.call(rbind, bands)
dat$px_or <- px(dat$OR); dat$px_lo <- px(pmax(orlo, dat$lo)); dat$px_hi <- px(pmin(orhi, dat$hi))
dat$msz <- scales::rescale(dat$wt, to = c(2.0, 3.5))
hy <- 1.7                                        # column-header row
ytop <- 2.7; ybot <- min(dat$y) - 0.7            # data extent
axy <- ybot - 1.0; favy <- ybot - 2.1            # axis + favors rows
brk <- c(0.5, 0.75, 1, 1.5, 2, 3)

f1 <- ggplot() +
  ## heavy top rule + header
  annotate("segment", x = 0, xend = 100, y = ytop, yend = ytop, colour = INKF, linewidth = 1.1) +
  annotate("text", x = 0,  y = hy, label = "Predictor",   hjust = 0, fontface = "bold", size = 3.1, colour = INKF, family = JBASE) +
  annotate("text", x = 67, y = hy, label = "OR (95% CI)", hjust = 0, fontface = "bold", size = 3.1, colour = INKF, family = JBASE) +
  annotate("text", x = 100, y = hy, label = "FDR", hjust = 1, fontface = "bold", size = 3.1, colour = INKF, family = JBASE) +
  annotate("segment", x = 0, xend = 100, y = hy - 0.75, yend = hy - 0.75, colour = INKF, linewidth = 0.5) +
  ## grey row-group bands + labels
  geom_rect(data = bnd, aes(xmin = 0, xmax = 100, ymin = y - 0.55, ymax = y + 0.55), fill = "#ECEFF1") +
  geom_text(data = bnd, aes(0, y, label = outcome), hjust = 0, fontface = "bold", size = 3.0, colour = INKF, family = JBASE) +
  ## reference line at OR = 1 (spans data rows only, not the header)
  annotate("segment", x = px(1), xend = px(1), y = ybot + 0.3, yend = hy - 0.95, colour = "#9AA0A6", linewidth = 0.5) +
  ## CI whiskers with end caps + square markers
  geom_segment(data = dat, aes(x = px_lo, xend = px_hi, y = y, yend = y), colour = INKF, linewidth = 0.5) +
  geom_segment(data = dat, aes(x = px_lo, xend = px_lo, y = y - 0.22, yend = y + 0.22), colour = INKF, linewidth = 0.5) +
  geom_segment(data = dat, aes(x = px_hi, xend = px_hi, y = y - 0.22, yend = y + 0.22), colour = INKF, linewidth = 0.5) +
  geom_point(data = dat, aes(px_or, y, size = msz), shape = 15, colour = INKF) +
  scale_size_identity() +
  ## text columns
  geom_text(data = dat, aes(3, y, label = index), hjust = 0, size = 2.9, colour = INKF, family = JBASE) +
  geom_text(data = dat, aes(67, y, label = ci), hjust = 0, size = 2.9, colour = INKF, family = JBASE,
            fontface = ifelse(dat$fdr_sig, "bold", "plain")) +
  geom_text(data = dat, aes(100, y, label = fdrv, fontface = ifelse(fdr_sig, "bold", "plain")),
            hjust = 1, size = 2.9, colour = INKF, family = JBASE) +
  ## heavy bottom rule
  annotate("segment", x = 0, xend = 100, y = ybot + 0.1, yend = ybot + 0.1, colour = INKF, linewidth = 1.1) +
  ## x-axis ticks + labels + favors
  annotate("segment", x = px(brk), xend = px(brk), y = axy + 0.35, yend = axy + 0.7, colour = INKF, linewidth = 0.4) +
  annotate("text", x = px(brk), y = axy, label = brk, size = 2.7, colour = INKF, family = JBASE) +
  annotate("segment", x = px(0.62), xend = px(0.5), y = favy, yend = favy, colour = INKF,
           linewidth = 0.4, arrow = arrow(length = unit(0.1, "cm"))) +
  annotate("segment", x = px(1.6), xend = px(2), y = favy, yend = favy, colour = INKF,
           linewidth = 0.4, arrow = arrow(length = unit(0.1, "cm"))) +
  annotate("text", x = px(0.7),  y = favy - 0.55, label = "Lower severity",  hjust = 0.5, size = 2.5, colour = MUTE, family = JBASE) +
  annotate("text", x = px(1.75), y = favy - 0.55, label = "Higher severity", hjust = 0.5, size = 2.5, colour = MUTE, family = JBASE) +
  coord_cartesian(xlim = c(0, 100), ylim = c(favy - 1.1, ytop + 0.4), expand = FALSE, clip = "off") +
  labs(title = "Adjusted association of aneurysm morphometry with aSAH severity",
       subtitle = "Proportional-odds ordinal regression, multiple imputation (m = 20); one index per model, OR per 1 SD.") +
  theme_void(base_family = JBASE) +
  theme(plot.title = element_text(face = "bold", size = 12.5, colour = INK, margin = margin(b = 2)),
        plot.subtitle = element_text(size = 8.5, colour = MUTE, margin = margin(b = 8)),
        plot.caption = element_text(size = 6.8, colour = MUTE, hjust = 0),
        plot.caption.position = "plot", plot.title.position = "plot",
        plot.margin = margin(8, 8, 6, 8), plot.background = element_rect(fill = "white", colour = NA)) +
  labs(caption = "Each estimate is a separate single-index model adjusted for age, sex, smoking, hypertension, diabetes, statin use, and aneurysm circulation. Marker size ∝ precision.\nBold, FDR < 0.05. Complete-case estimates in Table 2. FDR, Benjamini–Hochberg false-discovery rate; SD, standard deviation.")
save_fig(f1, "../figures/Figure2_forest", width = 8.8, height = 7.6)
dir.create("../figures/individual", showWarnings = FALSE)
save_fig(f1, "../figures/individual/Figure2", width = 8.8, height = 7.6)

## ================= FIGURE 2 : PCA size vs shape ==========================
## (A) PCA loading biplot (size axis PC1 vs shape axis PC2); (B) scree;
## (C) adjusted ordinal ORs of each component -> severity in the JAMA
##     table-forest format (monochrome squares, capped whiskers, OR column).
pc_in <- d[, c("neck","maxdiam","volume","aspect_ratio","hw_ratio","size_ratio","neck_parent")]
cc_pc <- complete.cases(pc_in)
pca <- prcomp(pc_in[cc_pc, ], scale. = TRUE)
ve <- summary(pca)$importance[2, 1:7]
## orient PC1 so larger = bigger; PC2 so positive = more elongated (aspect ratio)
if (pca$rotation["maxdiam","PC1"] < 0) { pca$rotation[,"PC1"] <- -pca$rotation[,"PC1"]; pca$x[,"PC1"] <- -pca$x[,"PC1"] }
if (pca$rotation["aspect_ratio","PC2"] < 0) { pca$rotation[,"PC2"] <- -pca$rotation[,"PC2"]; pca$x[,"PC2"] <- -pca$x[,"PC2"] }
load <- as.data.frame(pca$rotation[, 1:2])
load$var <- c("Neck width","Max diameter","Volume","Aspect ratio",
              "Height/width","Size ratio","Neck/parent")
load$axis <- ifelse(abs(load$PC1) >= abs(load$PC2), "Size axis (PC1)", "Shape axis (PC2)")
## two clearly distinct, colourblind-safe hues so size vs shape reads at a glance
SIZE_COL <- "#2C6FA6"   # blue  = size axis (PC1)
SHAPE_COL <- "#DF8F44"  # gold  = shape axis (PC2)
AXIS_COLS <- c("Size axis (PC1)" = SIZE_COL, "Shape axis (PC2)" = SHAPE_COL)
rng <- max(abs(c(load$PC1, load$PC2))) * 1.45

## (A) biplot — colour = axis on which each index loads more strongly
## (axis titles + caption carry the key; no in-plot legend to avoid clutter)
p2a <- ggplot(load, aes(PC1, PC2, colour = axis)) +
  geom_hline(yintercept = 0, colour = HAIR, linetype = 2) +
  geom_vline(xintercept = 0, colour = HAIR, linetype = 2) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.15,"cm"), type = "closed"), linewidth = 0.7) +
  geom_text_repel(aes(label = var), size = 2.8, family = JBASE, fontface = "bold",
                  min.segment.length = 0, box.padding = 0.5, point.padding = 0.3,
                  max.overlaps = Inf, seed = 7, segment.color = "grey75", segment.size = 0.3) +
  scale_colour_manual(values = AXIS_COLS, name = NULL) +
  ## colour key in the empty lower-left corner (no collision with arrows)
  annotate("point", x = -rng*0.26, y = -rng*0.58, shape = 15, size = 2.4, colour = SIZE_COL) +
  annotate("text",  x = -rng*0.18, y = -rng*0.58, label = "Size",  hjust = 0, size = 2.6, colour = SIZE_COL, fontface = "bold", family = JBASE) +
  annotate("point", x = -rng*0.26, y = -rng*0.74, shape = 15, size = 2.4, colour = SHAPE_COL) +
  annotate("text",  x = -rng*0.18, y = -rng*0.74, label = "Shape", hjust = 0, size = 2.6, colour = "#B5761F", fontface = "bold", family = JBASE) +
  coord_equal(xlim = c(-rng*0.28, rng*1.12), ylim = c(-rng, rng)) +
  labs(title = "Morphology splits into size and shape",
       subtitle = "The 7 indices cluster onto two axes",
       x = sprintf("PC1 — size axis (%.0f%%)", 100*ve[1]),
       y = sprintf("PC2 — shape axis (%.0f%%)", 100*ve[2])) +
  theme(legend.position = "none",
        plot.subtitle = element_text(size = 8, colour = MUTE))

## (B) scree
scree <- data.frame(PC = factor(paste0("PC",1:7), levels=paste0("PC",1:7)), ve = as.numeric(ve))
p2b <- ggplot(scree, aes(PC, ve)) +
  geom_col(fill = JAMA[["slate"]], width = 0.72) +
  geom_text(aes(label = sprintf("%.0f%%",100*ve)), vjust = -0.5, size = 2.6, colour = BODY, family = JBASE) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult=c(0,0.16))) +
  labs(title = "Variance explained", x = NULL, y = "Proportion of variance")

## (C) component -> severity, JAMA table-forest -----------------------------
dp <- d[cc_pc, ]; dp$PC1 <- pca$x[,"PC1"]; dp$PC2 <- pca$x[,"PC2"]
pcres <- data.frame()
for (on in names(OUTS)) for (pc in c("PC1","PC2")) {
  f <- as.formula(paste0(OUTS[[on]], " ~ ", pc, " + ", paste(COVARS, collapse=" + ")))
  r <- fit1(f, pc, dp); if (is.null(r)) next
  pcres <- rbind(pcres, data.frame(outcome=on, comp=pc, OR=exp(r["b"]),
    lo=exp(r["b"]-1.96*r["se"]), hi=exp(r["b"]+1.96*r["se"]),
    p=2*pnorm(abs(r["b"]/r["se"]), lower.tail=FALSE)))
}
comp_lab <- c(PC1 = "PC1 (size axis)", PC2 = "PC2 (shape axis)")
pfx0 <- 34; pfx1 <- 66; porlo <- 0.5; porhi <- 2
ppx <- function(or) pfx0 + (log10(or)-log10(porlo))/(log10(porhi)-log10(porlo))*(pfx1-pfx0)
prow <- list(); py <- 0; pband <- list()
for (o in names(OUTS)) {
  py <- py - 1.15
  pband[[length(pband)+1]] <- data.frame(y = py, outcome = names(which(OUTS==OUTS[[o]]))[1])
  pband[[length(pband)]]$outcome <- c(mFS="modified Fisher", HH="Hunt & Hess", WFNS="WFNS")[[o]]
  for (pc in c("PC1","PC2")) {
    r <- pcres[pcres$outcome==o & pcres$comp==pc, ]; py <- py - 1
    prow[[length(prow)+1]] <- data.frame(y = py, lab = comp_lab[[pc]],
      OR = r$OR, lo = r$lo, hi = r$hi, sig = (r$lo > 1 | r$hi < 1),
      ci = sprintf("%.2f (%.2f–%.2f)", r$OR, r$lo, r$hi))
  }
  py <- py - 0.5
}
pdat <- do.call(rbind, prow); pbnd <- do.call(rbind, pband)
pdat$px_or <- ppx(pdat$OR); pdat$px_lo <- ppx(pmax(porlo, pdat$lo)); pdat$px_hi <- ppx(pmin(porhi, pdat$hi))
phy <- 1.5; pytop <- 2.4; pybot <- min(pdat$y) - 0.7; paxy <- pybot - 0.9; pfavy <- pybot - 1.9
pbrk <- c(0.5, 0.7, 1, 1.5, 2)
p2c <- ggplot() +
  annotate("segment", x = 0, xend = 100, y = pytop, yend = pytop, colour = INKF, linewidth = 1.0) +
  annotate("text", x = 0,  y = phy, label = "Component",   hjust = 0, fontface = "bold", size = 3.0, colour = INKF, family = JBASE) +
  annotate("text", x = 70, y = phy, label = "OR (95% CI)", hjust = 0, fontface = "bold", size = 3.0, colour = INKF, family = JBASE) +
  annotate("segment", x = 0, xend = 100, y = phy - 0.72, yend = phy - 0.72, colour = INKF, linewidth = 0.5) +
  geom_rect(data = pbnd, aes(xmin = 0, xmax = 100, ymin = y - 0.55, ymax = y + 0.55), fill = "#ECEFF1") +
  geom_text(data = pbnd, aes(0, y, label = outcome), hjust = 0, fontface = "bold", size = 2.9, colour = INKF, family = JBASE) +
  annotate("segment", x = ppx(1), xend = ppx(1), y = pybot + 0.3, yend = phy - 0.9, colour = "#9AA0A6", linewidth = 0.5) +
  geom_segment(data = pdat, aes(x = px_lo, xend = px_hi, y = y, yend = y), colour = INKF, linewidth = 0.5) +
  geom_segment(data = pdat, aes(x = px_lo, xend = px_lo, y = y - 0.24, yend = y + 0.24), colour = INKF, linewidth = 0.5) +
  geom_segment(data = pdat, aes(x = px_hi, xend = px_hi, y = y - 0.24, yend = y + 0.24), colour = INKF, linewidth = 0.5) +
  geom_point(data = pdat, aes(px_or, y, colour = sig), shape = 15, size = 2.6) +
  scale_colour_manual(values = c(`TRUE` = unname(JAMA["brick"]), `FALSE` = INKF), guide = "none") +
  geom_text(data = pdat, aes(3, y, label = lab, fontface = ifelse(sig, "bold", "plain")), hjust = 0, size = 2.8, colour = INKF, family = JBASE) +
  geom_text(data = pdat, aes(70, y, label = ci, fontface = ifelse(sig, "bold", "plain")), hjust = 0, size = 2.8, colour = INKF, family = JBASE) +
  annotate("segment", x = 0, xend = 100, y = pybot + 0.1, yend = pybot + 0.1, colour = INKF, linewidth = 1.0) +
  annotate("segment", x = ppx(pbrk), xend = ppx(pbrk), y = paxy + 0.32, yend = paxy + 0.62, colour = INKF, linewidth = 0.4) +
  annotate("text", x = ppx(pbrk), y = paxy, label = pbrk, size = 2.6, colour = INKF, family = JBASE) +
  annotate("segment", x = ppx(0.62), xend = ppx(0.5), y = pfavy, yend = pfavy, colour = INKF, linewidth = 0.4, arrow = arrow(length = unit(0.09,"cm"))) +
  annotate("segment", x = ppx(1.6),  xend = ppx(2),   y = pfavy, yend = pfavy, colour = INKF, linewidth = 0.4, arrow = arrow(length = unit(0.09,"cm"))) +
  annotate("text", x = ppx(0.71), y = pfavy - 0.55, label = "Lower severity",  hjust = 0.5, size = 2.4, colour = MUTE, family = JBASE) +
  annotate("text", x = ppx(1.72), y = pfavy - 0.55, label = "Higher severity", hjust = 0.5, size = 2.4, colour = MUTE, family = JBASE) +
  annotate("text", x = 0, y = pytop + 1.0, hjust = 0, vjust = 1,
           label = "Size-axis (PC1) scores track severity most clearly for WFNS; shape-axis (PC2) scores do not (red = 95% CI excludes 1)",
           size = 2.7, colour = MUTE, family = JBASE) +
  coord_cartesian(xlim = c(0, 100), ylim = c(pfavy - 1.0, pytop + 2.0), expand = FALSE, clip = "off") +
  labs(title = "Which axis tracks severity?") +
  theme_void(base_family = JBASE) +
  theme(plot.title = element_text(face = "bold", size = 10.5, colour = INK, margin = margin(b = 2)),
        plot.title.position = "plot")

f2 <- (p2a | p2b) / p2c +
  plot_layout(heights = c(1.25, 1)) +
  plot_annotation(tag_levels = "A",
    caption = "A, PCA loading biplot; arrow colour marks the axis (blue = size/PC1, gold = shape/PC2) on which each index loads more strongly.\nC, Proportional-odds OR per 1 SD, adjusted for the same seven covariates as Figure 2. Marker and whisker style as in Figure 2; red = 95% CI excludes 1.",
    theme = theme_jama(base_size = 9) +
      theme(plot.caption = element_text(size = 6.8, colour = MUTE, hjust = 0))) &
  theme(plot.tag = element_text(face = "bold", size = 12, colour = INK))
save_fig(f2, "../figures/Figure3_PCA", width = 8.8, height = 8.4)
save_fig(p2a, "../figures/individual/Figure3a", width = 4.7, height = 4.7)
save_fig(p2b, "../figures/individual/Figure3b", width = 4.5, height = 4.2)
save_fig(p2c, "../figures/individual/Figure3c", width = 8.6, height = 4.4)

## ================= FIGURE 3 : sensitivity ================================
## (A) E-values for the FIVE FDR-significant associations (so every point shown
##     is genuinely FDR-significant; matches Table 3).
prim <- fp %>% filter(method=="Multiple imputation", p_FDR < 0.05)
ev <- do.call(rbind, lapply(seq_len(nrow(prim)), function(i){
  e <- EValue::evalues.OR(prim$OR[i], prim$lo[i], prim$hi[i], rare = FALSE)
  data.frame(label = paste(prim$outcome[i], prim$index[i], sep=" · "),
             point = e["E-values","point"], ci = e["E-values","lower"])
}))
ev$ci[is.na(ev$ci)] <- 1
ev$label <- reorder(ev$label, ev$ci)
p3a <- ggplot(ev, aes(ci, label)) +
  geom_vline(xintercept = 1, linetype = 2, color = MUTE, linewidth = 0.4) +
  geom_segment(aes(x = 1, xend = ci, y = label, yend = label), colour = HAIR, linewidth = 1.3) +
  geom_point(colour = unname(JAMA["brick"]), size = 3.1) +
  geom_text(aes(label = sprintf("%.2f", ci)), hjust = -0.45, size = 2.6, colour = INK, family = JBASE) +
  scale_x_continuous(limits = c(1, 1.62), breaks = seq(1, 1.6, 0.1),
                     expand = expansion(mult = c(0.01, 0.02))) +
  labs(title = "Robustness to unmeasured confounding (E-value)",
       subtitle = "How strong an unmeasured confounder — associated (risk ratio) with both morphology and severity —\nwould need to be to render each result non-significant. Farther from 1 = more robust.",
       x = "E-value  (95% CI limit, risk-ratio scale)", y = NULL) +
  theme(plot.subtitle = element_text(size = 7.3, colour = MUTE))

## (B) size-ratio tertiles -> HH probability of high grade
d$SRq <- cut(d$size_ratio, breaks = quantile(d$size_ratio, c(0,.33,.67,1), na.rm=TRUE),
             labels = c("Low","Mid","High"), include.lowest = TRUE)
thr <- d %>% filter(!is.na(SRq), !is.na(HH)) %>% group_by(SRq) %>%
  summarise(p_hi = mean(HH>=3), n = n(),
            se = sqrt(p_hi*(1-p_hi)/n), .groups="drop")
p3b <- ggplot(thr, aes(SRq, p_hi)) +
  geom_col(fill = "#AEBEC5", width = 0.6) +
  geom_errorbar(aes(ymin = pmax(0,p_hi-1.96*se), ymax = pmin(1,p_hi+1.96*se)), width = 0.16, colour = INK, linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f%% (n=%d)",100*p_hi,n)), vjust = -0.9, size = 2.6, colour = BODY, family = JBASE) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1), expand = expansion(mult=c(0,0.06))) +
  labs(title = "High-grade Hunt & Hess by size-ratio tertile", x = "Size-ratio tertile",
       y = "High-grade Hunt & Hess (3–5)")

## (C) MNAR delta robustness for size-ratio -> HH
deltas <- c(0, 0.5, 1, 1.5)
mn <- data.frame()
for (dl in deltas) {
  ests<-c(); vars<-c()
  for (k in 1:imp$m) {
    di <- complete(imp,k); miss <- is.na(d$HH_ord)
    hh <- as.integer(di$HH_ord); hh[miss] <- pmin(5, hh[miss]+dl)
    di$HHc <- factor(round(hh), ordered=TRUE)
    r <- fit1(as.formula(paste0("HHc ~ z_size_ratio + ", paste(COVARS, collapse=" + "))), "z_size_ratio", di)
    if(!is.null(r)){ ests<-c(ests,r["b"]); vars<-c(vars,r["se"]^2) }
  }
  m<-length(ests); qbar<-mean(ests); se<-sqrt(mean(vars)+(1+1/m)*var(ests))
  mn <- rbind(mn, data.frame(delta=dl, OR=exp(qbar), lo=exp(qbar-1.96*se), hi=exp(qbar+1.96*se)))
}
prim_or <- mn$OR[mn$delta == 0]
p3c <- ggplot(mn, aes(delta, OR)) +
  ## band = the primary (as-imputed) estimate; all shifts stay inside it
  geom_hline(yintercept = prim_or, colour = "#AEBEC5", linewidth = 3.2, alpha = 0.55) +
  geom_hline(yintercept = 1, linetype = 2, color = MUTE, linewidth = 0.4) +
  annotate("text", x = 1.5, y = 1.0, label = "no association", hjust = 1, vjust = -0.4, size = 2.1, colour = MUTE, family = JBASE) +
  annotate("text", x = 0, y = prim_or, label = "primary estimate", hjust = 0, vjust = -1.2, size = 2.1, colour = "#5B6B73", family = JBASE) +
  geom_line(colour = JAMA[["slate"]], linewidth = 0.5) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.06, color = JAMA[["slate"]], linewidth = 0.5) +
  geom_point(size = 2.5, shape = 15, color = JAMA[["slate"]]) +
  scale_y_continuous(limits = c(0.95, NA), breaks = seq(1, 2, 0.5),
                     expand = expansion(mult = c(0.02, 0.08))) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5),
                     labels = c("0\nas imputed", "0.5", "1.0", "1.5")) +
  labs(title = "If patients with missing grades were sicker",
       subtitle = "The size-ratio association stays essentially unchanged",
       x = "Assumed extra severity of patients with missing grades (grades)",
       y = "Adjusted OR per 1 SD") +
  theme(plot.margin = margin(6, 12, 6, 6),
        plot.subtitle = element_text(size = 7.3, colour = MUTE))

f3 <- p3a / (p3b | p3c) + plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "A",
    caption = "A, E-value (95% CI limit) per FDR-significant association; for scale, the smoking–lung-cancer E-value is ≈ 10.\nB, Complete-case proportions for size ratio (n = 133); error bars, 95% CI. C, Missing-not-at-random sensitivity; multiple-imputation estimate, error bars 95% CI.",
    theme = theme_jama(base_size = 9) +
      theme(plot.caption = element_text(size = 6.8, colour = MUTE, hjust = 0))) &
  theme(plot.tag = element_text(face = "bold", size = 12, colour = INK))
save_fig(f3, "../figures/Figure4_sensitivity", width = 9, height = 8.2)
save_fig(p3a, "../figures/individual/Figure4a", width = 8.4, height = 3.8)
save_fig(p3b, "../figures/individual/Figure4b", width = 4.4, height = 4.2)
save_fig(p3c, "../figures/individual/Figure4c", width = 4.7, height = 4.2)

message("Wrote Figure1_forest, Figure2_PCA, Figure3_sensitivity (.tiff + .pdf + .png)")
