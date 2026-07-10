## _journal_theme.R
## Shared JAMA house-style design system for all manuscript figures.
## Sourced by 04_figures.R and 06_method_figure.R so every panel shares one
## palette, one typeface, and one set of spacing rules.
##
## Palette: JAMA (ggsci::pal_jama) — the journal's own figure palette, chosen
## over the bright Okabe-Ito screen palette used previously.

suppressMessages({ library(ggplot2); library(grid) })

## ---- typeface ------------------------------------------------------------
## Helvetica Neue is JAMA/NEJM figure house type; fall back gracefully.
.pick_font <- function() {
  fam <- tryCatch(systemfonts::system_fonts()$family, error = function(e) character())
  for (f in c("Helvetica Neue", "Helvetica", "Arial")) if (f %in% fam) return(f)
  ""
}
JBASE <- .pick_font()

## ---- JAMA palette --------------------------------------------------------
JAMA <- c(
  slate  = "#374E55",  # primary ink / structural
  gold   = "#DF8F44",  # accent
  cyan   = "#00A1D5",  # nominal significance
  brick  = "#B24745",  # FDR significance / emphasis
  sage   = "#79AF97",  # positive / secondary series
  purple = "#6A6599",  # tertiary series
  taupe  = "#80796B"   # muted
)
## structural greys
INK    <- "#1A2A33"    # near-black title ink (warm slate)
BODY   <- "#374E55"    # body text = JAMA slate
MUTE   <- "#5B6B73"    # secondary text
HAIR   <- "#D4DBDE"    # hairline rules / gridlines
PANEL  <- "#FFFFFF"
NSGREY <- "#B7C0C4"    # non-significant markers

## significance colour map — restrained: slate + a single brick accent + grey.
## (Matches the austere methods figure; cyan/gold dropped as "un-JAMA".)
SIG_COLS <- c("FDR < 0.05" = unname(JAMA["brick"]),   # the one accent
              "p < 0.05"   = "#5B6B73",               # medium slate-grey
              "NS"         = NSGREY)                    # light grey

## ---- base theme ----------------------------------------------------------
theme_jama <- function(base_size = 9) {
  theme_minimal(base_size = base_size, base_family = JBASE) +
    theme(
      text             = element_text(colour = BODY),
      plot.title       = element_text(face = "bold", size = base_size + 2,
                                      colour = INK, margin = margin(b = 2)),
      plot.subtitle    = element_text(size = base_size - 0.5, colour = MUTE,
                                      margin = margin(b = 6)),
      plot.title.position = "plot",
      plot.caption     = element_text(size = base_size - 2, colour = MUTE, hjust = 0),
      plot.caption.position = "plot",
      axis.title       = element_text(size = base_size - 0.5, colour = MUTE),
      axis.text        = element_text(size = base_size - 1, colour = BODY),
      axis.ticks       = element_line(colour = HAIR, linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      panel.grid.major = element_line(colour = HAIR, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = PANEL, colour = NA),
      plot.background  = element_rect(fill = PANEL, colour = NA),
      strip.text       = element_text(face = "bold", size = base_size - 0.5, colour = INK),
      legend.position  = "bottom",
      legend.key.size  = unit(9, "pt"),
      legend.text      = element_text(size = base_size - 1, colour = BODY),
      legend.title     = element_text(size = base_size - 1, colour = MUTE),
      legend.margin    = margin(t = 2)
    )
}

## bold panel tag (A/B/C) styling for patchwork
jama_tag <- function() theme(
  plot.tag = element_text(face = "bold", size = 11, colour = INK),
  plot.tag.position = c(0, 1)
)

## device-safe save: TIFF (300 dpi, LZW), vector PDF, and ragg PNG preview.
save_fig <- function(plot, stem, width, height) {
  ragg_ok <- requireNamespace("ragg", quietly = TRUE)
  ggsave(paste0(stem, ".tiff"), plot, width = width, height = height,
         dpi = 300, compression = "lzw",
         device = if (ragg_ok) ragg::agg_tiff else grDevices::tiff, bg = "white")
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height,
         device = grDevices::cairo_pdf, bg = "white")
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 200,
         device = if (ragg_ok) ragg::agg_png else grDevices::png, bg = "white")
  invisible(stem)
}

theme_set(theme_jama())
