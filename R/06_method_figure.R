## 06_method_figure.R
## Analytic-workflow schematic (manuscript Figure 4 + repo). Pure ggplot ->
## reproducible, editable vector output. Colorblind-safe (Okabe-Ito).

suppressMessages({ library(ggplot2) })
dir.create("../figures", showWarnings = FALSE)

OK <- c(blue="#0072B2", orange="#E69F00", green="#009E73", red="#D55E00",
        purple="#CC79A7", sky="#56B4E9", grey="#7F7F7F", light="#EAF3FA")

boxes <- data.frame(
  id   = c("cohort","morph","out","miss","model","pca","ridge","logit",
           "mult","sens","result"),
  x    = c(5,    2.6,  7.4,  5,    5,    1.6,  1.6,  8.4,  8.4,  5,    5),
  y    = c(9.3,  7.8,  7.8,  6.3,  4.9,  4.9,  3.6,  4.9,  3.6,  2.3,  0.7),
  w    = c(6.4,  3.4,  3.4,  6.4,  3.6,  2.7,  2.7,  2.7,  2.7,  8.6,  7.2),
  h    = c(0.9,  1.15, 1.15, 0.95, 1.0,  1.0,  0.9,  1.0,  0.9,  1.15, 1.0),
  fill = c(OK["sky"], OK["light"], OK["light"], OK["orange"], OK["green"],
           OK["light"], OK["light"], OK["light"], OK["light"], OK["purple"], OK["blue"]),
  col  = c("white","black","black","black","white","black","black","black",
           "black","white","white"),
  lab = c(
    "aSAH cohort — 171 patients\nSingle center, 2013–2023",
    "7 morphometric indices\n(admission imaging)\nsize ratio · aspect ratio · max diameter\nvolume · neck · H/W · neck/parent",
    "Severity outcomes\nHunt & Hess · WFNS · mFS (ordinal)\nGCS (continuous)",
    "Missing data ~22% (grades, HTN)  →  Multiple imputation (MICE, m=20, primary)\n+ complete-case (sensitivity)",
    "PRIMARY MODEL\nProportional-odds ordinal\nregression, per-SD,\nadjusted (7 covariates)*",
    "PCA composite\nPC1 size / PC2 shape",
    "Ridge\n(ranking only)",
    "Logistic\n(dichotomized)",
    "BH-FDR\nmultiplicity",
    "Sensitivity suite:  PO assumption (Brant) · MI-then-delete · MNAR delta-shift · E-values · bootstrap AUC",
    "RESULT — size-related indices (size ratio, aspect ratio) predict higher aSAH severity;\nsignal loads on the size axis (PC1); robust to imputation & unmeasured confounding"),
  stringsAsFactors = FALSE
)

seg <- data.frame(
  x    = c(5,    5,    5,    5,    3.3,  6.7,  8.4,  5,    5),
  xend = c(3.3,  6.7,  5,    5,    1.6,  8.4,  8.4,  5,    5),
  y    = c(8.85, 8.85, 7.22, 5.82, 4.4,  4.4,  4.4,  4.4,  1.72),
  yend = c(8.38, 8.38, 6.78, 5.4,  4.1,  4.1,  4.1,  2.88, 1.2)
)

adjust_note <- "*All adjusted models control for age, sex, smoking, hypertension, diabetes, statin use, and aneurysm circulation."

p <- ggplot() +
  geom_segment(data = seg, aes(x, y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
               linewidth = 0.5, color = OK["grey"]) +
  geom_tile(data = boxes, aes(x, y, width = w, height = h, fill = I(fill)),
            color = "grey40", linewidth = 0.4) +
  geom_text(data = boxes, aes(x, y, label = lab, color = I(col)),
            size = 2.75, lineheight = 0.95, fontface = "plain") +
  annotate("text", x = 0.15, y = -0.25, label = adjust_note, size = 2.2,
           hjust = 0, fontface = "italic", color = OK["grey"]) +
  annotate("text", x = 5, y = 9.95, label = "Analytic workflow",
           size = 5, fontface = "bold") +
  coord_cartesian(xlim = c(0, 10), ylim = c(-0.6, 10.3), expand = FALSE) +
  theme_void() + theme(plot.margin = margin(6,6,6,6))

ggsave("../figures/Figure4_method.tiff", p, width = 9, height = 8.2, dpi = 300, compression = "lzw", type = "cairo")
ggsave("../figures/Figure4_method.pdf",  p, width = 9, height = 8.2)
message("Wrote Figure4_method (.tiff + .pdf)")
