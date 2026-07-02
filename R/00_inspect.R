## 00_inspect.R -- one-off: report coding of key variables so cleaning is correct.
## Not part of the reproducible pipeline; diagnostic only.
suppressMessages(library(readxl))

xlsx <- Sys.getenv("ASAH_XLSX",
  unset = "../aSAH_Comprehensive_Data_1.1.26.xlsx")

d <- suppressMessages(read_excel(xlsx, sheet = 1, skip = 1, .name_repair = "minimal"))
d <- as.data.frame(d)

show <- function(i, lab) {
  v <- d[[i]]
  tab <- sort(table(v, useNA = "ifany"), decreasing = TRUE)
  cat(sprintf("[%d] %-12s class=%-9s uniq=%d\n", i, lab, class(v)[1], length(tab)))
  cat("     top:", paste(head(names(tab), 12), collapse = " | "), "\n")
}

vars <- list(
  c(7, "Sex"), c(9, "Smoking"), c(10, "Diabetes"), c(11, "HTN"),
  c(13, "Statin"), c(23, "Circulation"), c(26, "HH"), c(27, "Fisher"),
  c(28, "WFNS"), c(29, "GCS"), c(40, "Bleb"), c(19, "Morphology"),
  c(37, "Neck"), c(38, "MaxDiam"), c(39, "Volume"), c(42, "AspectR"),
  c(43, "HWR"), c(44, "SizeR"), c(45, "NeckParent"), c(41, "ParentDiam")
)
for (x in vars) show(as.integer(x[[1]]), x[[2]])
