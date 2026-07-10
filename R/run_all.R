## run_all.R -- reproduce the entire analysis from the raw workbook.
## Usage:  ASAH_XLSX="/path/to/aSAH_Comprehensive_Data_1.1.26.xlsx" Rscript run_all.R
## Each step writes to ../outputs and ../figures.

stopifnot(nzchar(Sys.getenv("ASAH_XLSX")))
steps <- c("01_clean.R", "02_analysis.R", "03_rwe_checks.R",
           "05_tables.R", "07_tables_gt.R", "04_figures.R", "06_method_figure.R")
for (s in steps) {
  message("\n==== running ", s, " ====")
  source(s, local = new.env())
}
message("\nAll steps complete. See ../outputs/ and ../figures/.")
