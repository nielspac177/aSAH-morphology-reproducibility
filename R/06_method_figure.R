## 06_method_figure.R  — Figure 4: Analytic workflow (JAMA house style)
## Swim-lane methods schematic: rotated stage labels down the left margin,
## content boxes connected by arrows (across each lane, and down between lanes).
## Methods only. Reproducible vector + 300-dpi raster; palette from _journal_theme.R.

suppressMessages({ library(ggplot2); library(grid); library(ggtext) })
source("_journal_theme.R")
dir.create("../figures", showWarnings = FALSE)

BORD <- JAMA[["slate"]]     # box borders / arrows
ACC  <- JAMA[["brick"]]     # primary-model emphasis
LANEF <- "#EEF3F5"          # lane-label fill
BOXF <- "white"

## ---- helpers -------------------------------------------------------------
## sharp-cornered content box (grid rect so borders stay crisp at any aspect)
rectg <- function(xmin, xmax, ymin, ymax, fill = BOXF, col = BORD, lwd = 0.8)
  annotation_custom(rectGrob(gp = gpar(fill = fill, col = col, lwd = lwd)),
                    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
## rounded lane-label box
laneg <- function(ymin, ymax, xmin = 1.5, xmax = 8)
  annotation_custom(roundrectGrob(r = unit(4, "pt"),
      gp = gpar(fill = LANEF, col = BORD, lwd = 0.8)),
      xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
arr <- function() arrow(length = unit(0.13, "cm"), type = "closed")

## ---- lane definitions (top -> bottom) ------------------------------------
lanes <- data.frame(
  label = c("Cohort &\nmeasurements", "Data\npreparation", "Primary\nmodel",
            "Multiplicity", "Sensitivity &\nvalidation"),
  ymin  = c(84, 70, 46, 30, 14),
  ymax  = c(96, 82, 68, 44, 28))
lanes$yc <- (lanes$ymin + lanes$ymax) / 2

## ---- content boxes -------------------------------------------------------
boxes <- data.frame(
  xmin = c(11, 40, 71,        # L1
           11, 43, 78,        # L2
           11,                # L3 (wide, drawn separately)
           11, 40, 69,        # L4
           11, 34, 57, 80),   # L5
  xmax = c(31, 62, 98,
           34, 70, 98,
           98,
           37, 66, 98,
           31, 54, 77, 98),
  ymin = c(85.5,85.5,85.5,  71.5,71.5,71.5,  47.5,  31.5,31.5,31.5,  15.5,15.5,15.5,15.5),
  ymax = c(94.5,94.5,94.5,  80.5,80.5,80.5,  66.5,  42.5,42.5,42.5,  26.5,26.5,26.5,26.5),
  label = c(
    "171 patients\nsingle-center aSAH\n2013–2023",
    "7 morphometric indices\nadmission CTA / DSA",
    "4 severity outcomes\nHunt–Hess · WFNS · mFS (ordinal)\nGCS (continuous)",
    "Standardize predictors\n(per SD) · log-transform volume",
    "Multiple imputation\nMICE, m = 20 (~22% missing)",
    "Complete-case\nanalysis (sensitivity)",
    "",  # primary-model box filled by rich text below
    "Benjamini–Hochberg\nFDR control",
    "Brant test\n(parallel-odds)",
    "MI-then-delete",
    "E-values", "MNAR δ-shift", "Bootstrap-\nvalidated AUC", "PCA: size vs\nshape axis"),
  stringsAsFactors = FALSE)
prim <- 7L                                       # index of the wide primary box
boxes$xc <- (boxes$xmin + boxes$xmax)/2
boxes$yc <- (boxes$ymin + boxes$ymax)/2

p <- ggplot() +
  ## lane labels
  lapply(seq_len(nrow(lanes)), function(i) laneg(lanes$ymin[i], lanes$ymax[i])) +
  geom_text(data = lanes, aes(4.75, yc, label = label), angle = 90,
            fontface = "bold", size = 2.75, colour = BORD, lineheight = 0.95, family = JBASE) +

  ## content boxes (all except the wide primary box)
  lapply(setdiff(seq_len(nrow(boxes)), prim), function(i)
    rectg(boxes$xmin[i], boxes$xmax[i], boxes$ymin[i], boxes$ymax[i])) +
  ## primary-model box (emphasized border)
  rectg(boxes$xmin[prim], boxes$xmax[prim], boxes$ymin[prim], boxes$ymax[prim],
        fill = "#FBF6F5", col = ACC, lwd = 1.5) +

  ## box text (multi-line, centered)
  geom_text(data = boxes[-prim, ], aes(xc, yc, label = label),
            size = 2.5, colour = BODY, lineheight = 1.05, family = JBASE) +

  ## ---- primary-model rich content ----
  annotate("text", x = 13, y = 64.7, label = "Primary model — proportional-odds ordinal regression",
           hjust = 0, fontface = "bold", size = 3.05, colour = ACC, family = JBASE) +
  annotate("segment", x = 13, xend = 96, y = 62.3, yend = 62.3, colour = "#E2C9C7", linewidth = 0.4) +
  geom_richtext(aes(54.5, 58.6,
      label = "logit&#8202;[&#8202;Pr(&#8202;*Y*&#8202;≥&#8202;*k*&#8202;)&#8202;] = α<sub>*k*</sub> + β·*Z* + Σ<sub>*j*</sub> γ<sub>*j*</sub>&#8202;*C<sub>j</sub>*"),
      hjust = 0.5, size = 3.4, colour = INK, fill = NA, label.color = NA, family = JBASE) +
  geom_richtext(aes(13, 54,
      label = "*Z* = morphometric index (per 1 SD) · OR per 1 SD = *e*<sup>β</sup> · OR > 1 = higher severity · one index per model"),
      hjust = 0, size = 2.5, colour = BODY, fill = NA, label.color = NA, family = JBASE) +
  geom_richtext(aes(13, 50.4,
      label = "<span style='color:#5B6B73'>Adjusted for</span> **age · sex · smoking · hypertension · diabetes · statin use · aneurysm circulation**"),
      hjust = 0, size = 2.5, colour = BODY, fill = NA, label.color = NA, family = JBASE) +

  ## ---- horizontal arrows (within lanes) ----
  annotate("segment", x = c(31,62, 34,70, 31,54,77, 37,66),
                     xend = c(40,71, 43,78, 34,57,80, 40,69),
                     y   = c(90,90, 76,76, 21,21,21, 37,37),
                     yend= c(90,90, 76,76, 21,21,21, 37,37),
           colour = BORD, linewidth = 0.6, arrow = arr()) +
  ## ---- vertical arrows (between lanes), down the left column ----
  annotate("segment", x = 21, xend = 21,
           y    = c(85.5, 71.5, 47.5, 31.5),
           yend = c(80.5, 66.7, 42.5, 26.5),
           colour = BORD, linewidth = 0.6, arrow = arr()) +

  ## ---- title + footnote ----
  annotate("text", x = 1.5, y = 99.4, label = "Analytic workflow", hjust = 0,
           fontface = "bold", size = 4.6, colour = INK, family = JBASE) +
  annotate("segment", x = 1.5, xend = 98, y = 8.5, yend = 8.5, colour = HAIR, linewidth = 0.35) +
  geom_textbox(aes(1.5, 6.5,
      label = "Single-center retrospective cohort; each morphometric index modelled separately per outcome. Continuous outcome (admission GCS; lower = worse) modelled by linear regression. CTA, CT angiography; DSA, digital subtraction angiography; WFNS, World Federation of Neurosurgical Societies grade; mFS, modified Fisher scale; FDR, false-discovery rate; MNAR, missing not at random; SD, standard deviation."),
      hjust = 0, vjust = 1, halign = 0, width = unit(0.97, "npc"),
      size = 1.95, colour = MUTE, fill = NA, box.color = NA,
      box.padding = unit(c(0,0,0,0),"pt"), lineheight = 1.18, family = JBASE) +

  coord_cartesian(xlim = c(0, 100), ylim = c(2, 100.5), expand = FALSE, clip = "off") +
  theme_void() +
  theme(plot.margin = margin(7, 9, 7, 9),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA))

save_fig(p, "../figures/Figure1_method", width = 8.7, height = 10.0)
## also export as a standalone (single-panel) tiff for per-figure submission
dir.create("../figures/individual", showWarnings = FALSE)
save_fig(p, "../figures/individual/Figure1", width = 8.7, height = 10.0)
message("Wrote Figure1_method (.tiff + .pdf + .png)")
