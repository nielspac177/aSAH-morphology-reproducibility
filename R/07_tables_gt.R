## 07_tables_gt.R  — Publication tables typeset with gt (JAMA house style).
## Reads the CSVs written by 05_tables.R and renders Tables 1–3 + Supp S2 as
## PNG (300-dpi-equiv), vector PDF, and HTML. Horizontal hairline rules only,
## Helvetica, bold column spanners, FDR-significant rows emphasized in brick.

suppressMessages({ library(gt); library(dplyr) })
td <- "../outputs/tables"; dir.create(td, showWarnings = FALSE, recursive = TRUE)

SLATE <- "#374E55"; BRICK <- "#B24745"; MUTE <- "#5B6B73"; HAIR <- "#D4DBDE"
FONT  <- c("Helvetica Neue", "Helvetica", "Arial", gt::default_fonts())

## shared JAMA styling ------------------------------------------------------
jama_gt <- function(g, note = NULL) {
  g <- g |>
    opt_table_font(font = FONT) |>
    tab_options(
      table.font.size = px(13),
      table.font.color = SLATE,
      heading.title.font.size = px(15),
      heading.title.font.weight = "bold",
      heading.subtitle.font.size = px(12),
      heading.align = "left",
      table.border.top.style = "solid", table.border.top.width = px(2), table.border.top.color = SLATE,
      table.border.bottom.style = "solid", table.border.bottom.width = px(2), table.border.bottom.color = SLATE,
      column_labels.border.top.style = "none",
      column_labels.border.bottom.style = "solid", column_labels.border.bottom.width = px(1.5), column_labels.border.bottom.color = SLATE,
      column_labels.font.weight = "bold",
      row_group.border.top.style = "solid", row_group.border.top.width = px(1), row_group.border.top.color = HAIR,
      row_group.border.bottom.style = "none",
      row_group.font.weight = "bold",
      table_body.border.bottom.color = SLATE,
      table_body.hlines.style = "none",
      table_body.vlines.style = "none",
      stub.border.style = "none",
      data_row.padding = px(4),
      column_labels.padding = px(6),
      source_notes.font.size = px(10.5), source_notes.padding = px(6)
    )
  if (!is.null(note)) g <- g |> tab_source_note(source_note = note)
  g
}
save_tab <- function(g, stem, vwidth = 1000) {
  gt::gtsave(g, paste0(stem, ".png"), expand = 8, zoom = 3, vwidth = vwidth)
  gt::gtsave(g, paste0(stem, ".html"))
  tryCatch(gt::gtsave(g, paste0(stem, ".pdf"), vwidth = vwidth),
           error = function(e) message("pdf skipped: ", conditionMessage(e)))
}

## ======================= TABLE 1 =========================================
t1 <- read.csv(file.path(td, "Table1_descriptives.csv"), check.names = FALSE)
t1 <- t1[t1$Variable != "N", ]
g1 <- t1 |>
  gt() |>
  tab_header(title = "Table 1. Cohort characteristics",
             subtitle = "Aneurysmal subarachnoid hemorrhage, N = 171") |>
  cols_label(Variable = "Characteristic", Value = "Value") |>
  cols_align("left", Variable) |>
  cols_align("right", Value) |>
  jama_gt(note = md("Values are mean (SD) or n (%). HH, Hunt & Hess; mFS, modified Fisher scale; WFNS, World Federation of Neurosurgical Societies grade. Percentages use non-missing denominators."))
save_tab(g1, file.path(td, "Table1"), vwidth = 620)

## ======================= TABLE 2 =========================================
t2 <- read.csv(file.path(td, "Table2_ordinal_ORs.csv"), check.names = FALSE)
t2$fdr_fmt <- ifelse(t2$mi_fdr < 0.001, "<0.001", sprintf("%.3f", t2$mi_fdr))
t2$sig <- t2$mi_fdr < 0.05
ord_out <- c("modified Fisher", "Hunt & Hess", "WFNS")
t2$outcome <- factor(t2$outcome, levels = ord_out)
t2 <- t2[order(t2$outcome), ]
g2 <- t2 |>
  select(outcome, index, N, cc, mi, fdr_fmt, sig) |>
  gt(groupname_col = "outcome", rowname_col = "index") |>
  tab_header(title = "Table 2. Adjusted per-SD ordinal odds ratios for morphology and aSAH severity",
             subtitle = "Proportional-odds ordinal regression, one index per model; OR per 1 SD, OR > 1 = higher severity") |>
  cols_label(N = "No.", cc = md("**Complete-case**<br>OR (95% CI), *p*"),
             mi = md("**Multiple imputation**<br>OR (95% CI), *p*"), fdr_fmt = "MI FDR") |>
  cols_align("center", c(N, fdr_fmt)) |>
  cols_align("left", c(cc, mi)) |>
  tab_stubhead(label = "Predictor") |>
  tab_style(style = cell_text(color = BRICK, weight = "bold"),
            locations = cells_body(columns = fdr_fmt, rows = sig)) |>
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_stub(rows = sig)) |>
  tab_style(style = list(cell_fill(color = "#ECEFF1"), cell_text(weight = "bold")),
            locations = cells_row_groups()) |>
  cols_hide(sig) |>
  jama_gt(note = md("Values are proportional-odds (cumulative) odds ratios per 1 SD from separate single-index models, adjusted for age, sex, smoking, hypertension, diabetes, statin use, and aneurysm circulation. Multiple imputation: MICE, m = 20. Benjamini–Hochberg FDR controls across the 21 ordinal tests (7 indices × 3 grades); **bold** = FDR < 0.05. The continuous outcome (admission GCS) was analysed separately by linear regression. SD, standard deviation."))
save_tab(g2, file.path(td, "Table2"), vwidth = 1040)

## ======================= TABLE 3 =========================================
t3 <- read.csv(file.path(td, "Table3_FDR_significant.csv"), check.names = FALSE)
t3$fdr_fmt <- ifelse(t3$mi_fdr < 0.001, "<0.001", sprintf("%.3f", t3$mi_fdr))
g3 <- t3 |>
  select(outcome, index, mi, fdr_fmt) |>
  gt() |>
  tab_header(title = "Table 3. Associations surviving FDR correction",
             subtitle = "Multiple-imputation estimates, ranked by FDR") |>
  cols_label(outcome = "Outcome", index = "Predictor",
             mi = md("OR (95% CI), *p*"), fdr_fmt = "FDR") |>
  cols_align("left", c(outcome, index, mi)) |>
  cols_align("center", fdr_fmt) |>
  tab_style(style = cell_text(color = BRICK, weight = "bold"),
            locations = cells_body(columns = fdr_fmt)) |>
  jama_gt(note = md("Proportional-odds (cumulative) odds ratios per 1 SD, adjusted, from multiple imputation. Five of the 21 ordinal morphology–severity tests (7 indices × 3 grades) survive Benjamini–Hochberg correction."))
save_tab(g3, file.path(td, "Table3"), vwidth = 720)

message("Wrote Table1, Table2, Table3 (.png + .html + .pdf) to ", td)
