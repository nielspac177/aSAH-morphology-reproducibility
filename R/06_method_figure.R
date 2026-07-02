## 06_method_figure.R
## Analytic-workflow schematic (manuscript Figure 4 + repo). Staged-pipeline
## design: numbered rail, rounded content cards, highlighted primary-model stage,
## result banner. Reproducible vector output; colorblind-safe, restrained palette.

suppressMessages({ library(ggplot2); library(grid) })
dir.create("../figures", showWarnings = FALSE)

navy   <- "#123B5E"; ink    <- "#16324F"; sub <- "#33455A"
boxf   <- "#EEF3F8"; boxc   <- "#B9C8D8"
hlf    <- "#E3F1EC"; hlc    <- "#4FA98A"
resf   <- "#123B5E"; foot   <- "#7A8896"

stages <- data.frame(
  n = 1:5,
  y = c(8.8, 7.3, 5.8, 4.3, 2.8),
  title = c("Cohort & measurements", "Data preparation", "Primary model",
            "Secondary models & multiplicity", "Robustness"),
  content = c(
    "171 patients · single center, 2013–2023\n7 morphometric indices (admission imaging)\n4 outcomes: Hunt–Hess · WFNS · mFS (ordinal) · GCS (continuous)",
    "standardize predictors (per SD) · log-transform volume\nmultiple imputation for ~22% missing (MICE, m = 20)\ncomplete-case analysis retained as sensitivity",
    "proportional-odds ordinal regression\none index per model, adjusted for 7 covariates*\neffect = odds ratio per 1 SD  (> 1 → higher severity)",
    "PCA: size axis (PC1) vs shape axis (PC2)\nridge (ranking only) · dichotomized logistic\nBenjamini–Hochberg false-discovery-rate control",
    "proportional-odds check (Brant) · MI-then-delete\nMNAR delta-shift · E-values · bootstrap-validated AUC"),
  hl = c(FALSE, FALSE, TRUE, FALSE, FALSE),
  stringsAsFactors = FALSE)

bx0 <- 1.75; bx1 <- 9.6; hh <- 0.62
rr <- function(ymid, fill, col, lwd = 1.2, xmin = bx0, xmax = bx1, half = hh)
  annotation_custom(roundrectGrob(r = unit(6, "pt"),
    gp = gpar(fill = fill, col = col, lwd = lwd)),
    xmin = xmin, xmax = xmax, ymin = ymid - half, ymax = ymid + half)

boxes <- lapply(seq_len(nrow(stages)), function(i)
  rr(stages$y[i], ifelse(stages$hl[i], hlf, boxf), ifelse(stages$hl[i], hlc, boxc),
     lwd = ifelse(stages$hl[i], 1.6, 1.2)))

res_y <- 1.2
banner <- rr(res_y, resf, NA, xmin = 0.4, xmax = 9.6, half = 0.56)

p <- ggplot() +
  # vertical rail behind the cards
  annotate("segment", x = 0.95, xend = 0.95, y = stages$y[1], yend = stages$y[5],
           color = navy, linewidth = 1.1) +
  boxes + banner +
  # rail continues into the result banner with an arrowhead
  annotate("segment", x = 0.95, xend = 0.95, y = stages$y[5] - hh, yend = res_y + 0.62,
           color = navy, linewidth = 1.1,
           arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  # numbered stage nodes
  geom_point(data = stages, aes(0.95, y), size = 9, shape = 21,
             fill = navy, color = "white", stroke = 1.1) +
  geom_text(data = stages, aes(0.95, y, label = n), color = "white",
            fontface = "bold", size = 3.7) +
  # card titles + content
  geom_text(data = stages, aes(bx0 + 0.18, y + 0.42, label = title), hjust = 0,
            fontface = "bold", size = 3.35, color = ink) +
  geom_text(data = stages, aes(bx0 + 0.18, y + 0.18, label = content), hjust = 0,
            vjust = 1, size = 2.62, color = sub, lineheight = 1.08) +
  # result banner text
  annotate("text", x = 5, y = res_y + 0.17,
           label = "Size ratio and aspect ratio — the size axis — track greater aSAH severity",
           color = "white", fontface = "bold", size = 3.15) +
  annotate("text", x = 5, y = res_y - 0.2,
           label = "consistent across grading scales  ·  robust to imputation and unmeasured confounding",
           color = "#D5E1EC", size = 2.72) +
  # title + footnote
  annotate("text", x = 0.4, y = 9.9, label = "Analytic workflow", hjust = 0,
           fontface = "bold", size = 5.4, color = ink) +
  annotate("text", x = 0.4, y = 0.25,
           label = "*Adjusted for age, sex, smoking, hypertension, diabetes, statin use, and aneurysm circulation.",
           hjust = 0, fontface = "italic", size = 2.4, color = foot) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 10.25), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(8, 10, 8, 10),
        plot.background  = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

ggsave("../figures/Figure4_method.tiff", p, width = 8.6, height = 9.4, dpi = 300, compression = "lzw", type = "cairo")
ggsave("../figures/Figure4_method.pdf",  p, width = 8.6, height = 9.4)
ggsave("../figures/Figure4_method.png",  p, width = 8.6, height = 9.4, dpi = 150, type = "cairo")
message("Wrote Figure4_method (.tiff + .pdf + .png)")
