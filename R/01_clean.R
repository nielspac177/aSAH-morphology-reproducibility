## 01_clean.R
## Build the analysis dataset from the raw clinical workbook.
## Reads the raw xlsx (path via ASAH_XLSX env var; kept OUT of the repo).
## Writes a de-identified analysis data frame to outputs/analysis_data.rds (git-ignored).
##
## Ground truth for coding decisions = manuscript Table 1. The tail of this
## script prints a reproduction of Table 1 so any recoding error is caught.

suppressMessages({
  library(readxl)
})

xlsx <- Sys.getenv("ASAH_XLSX")
if (!nzchar(xlsx)) stop("Set ASAH_XLSX to the raw workbook path.")

raw <- suppressWarnings(suppressMessages(
  read_excel(xlsx, sheet = 1, skip = 1, .name_repair = "minimal")))
raw <- as.data.frame(raw)

## --- helpers -------------------------------------------------------------
num  <- function(x) suppressWarnings(as.numeric(trimws(as.character(x))))
blank_na <- function(x) { x <- trimws(as.character(x)); x[x %in% c("", "NA", "N/A", "na")] <- NA; x }

## --- pull columns by position (header text is unstable) ------------------
d <- data.frame(
  id          = seq_len(nrow(raw)),
  sex_raw     = blank_na(raw[[7]]),
  age         = num(raw[[8]]),
  smoke_raw   = blank_na(raw[[9]]),
  dm_raw      = blank_na(raw[[10]]),
  htn_raw     = blank_na(raw[[11]]),
  statin_raw  = blank_na(raw[[13]]),
  circ_raw    = blank_na(raw[[23]]),
  location    = blank_na(raw[[22]]),
  HH          = num(raw[[26]]),
  mFS         = num(raw[[27]]),   # "Fisher Grade" col = modified Fisher in ms
  WFNS        = num(raw[[28]]),
  GCS         = num(raw[[29]]),
  bleb        = num(raw[[40]]),
  neck        = num(raw[[37]]),
  maxdiam     = num(raw[[38]]),
  volume      = num(raw[[39]]),
  aspect_ratio= num(raw[[42]]),   # dome height / neck width -- NOT "ASPECTS"
  hw_ratio    = num(raw[[43]]),
  size_ratio  = num(raw[[44]]),
  neck_parent = num(raw[[45]]),
  parent_diam = num(raw[[41]]),
  stringsAsFactors = FALSE
)

## --- harmonize covariates ------------------------------------------------
# Sex: workbook mixes {Female,Male} with numeric {0,1}. The numeric block is coded
# 0 = Female, 1 = Male (this reproduces the manuscript's 113 F / 58 M exactly).
d$sex <- with(d, ifelse(grepl("^f", sex_raw, ignore.case = TRUE) | sex_raw == "0",
                        "Female",
                 ifelse(grepl("^m", sex_raw, ignore.case = TRUE) | sex_raw == "1",
                        "Male", NA)))

# Smoking: collapse the 9 messy codings to Never / Former / Current (Table 1: 92/22/51).
smk <- tolower(d$smoke_raw)
d$smoking <- NA_character_
d$smoking[smk %in% c("0", "never smoker", "never tobacco user", "unknown if ever smoked", "never")] <- "Never"
d$smoking[smk %in% c("1", "former", "former smoker", "former tobacco user")] <- "Former"
d$smoking[grepl("^2", smk) | smk %in% c("current", "current smoker")] <- "Current"
# binary version used for adjustment (fewer empty cells -> avoids rank deficiency)
d$smoking_ever <- ifelse(is.na(d$smoking), NA, ifelse(d$smoking == "Never", "Never", "Ever"))

bin01 <- function(x) { x <- as.character(x); ifelse(x == "1", "Yes", ifelse(x == "0", "No", NA)) }
d$dm     <- bin01(d$dm_raw)
d$htn    <- bin01(d$htn_raw)
d$statin <- bin01(d$statin_raw)
d$circulation <- ifelse(d$circ_raw == "1", "Posterior",
                 ifelse(d$circ_raw == "0", "Anterior", NA))

## --- outcomes ------------------------------------------------------------
d$HH_ord   <- factor(d$HH,   levels = 1:5, ordered = TRUE)
d$mFS_ord  <- factor(d$mFS,  levels = 1:4, ordered = TRUE)
d$WFNS_ord <- factor(d$WFNS, levels = 1:5, ordered = TRUE)

d$HH_hi   <- ifelse(d$HH   >= 3, 1L, 0L)   # HH 3-5 vs 1-2
d$mFS_hi  <- ifelse(d$mFS  >= 3, 1L, 0L)   # mFS 3-4 vs 1-2
d$WFNS_hi <- ifelse(d$WFNS >= 4, 1L, 0L)   # WFNS 4-5 vs 1-3

## --- transformed / standardized predictors ------------------------------
d$log_volume <- log(d$volume + 1)          # volume is severely right-skewed
morph <- c("neck", "maxdiam", "volume", "log_volume",
           "aspect_ratio", "hw_ratio", "size_ratio", "neck_parent")
for (m in morph) d[[paste0("z_", m)]] <- as.numeric(scale(d[[m]]))

## --- factor covariates ---------------------------------------------------
d$sex          <- factor(d$sex)
d$smoking      <- factor(d$smoking, levels = c("Never", "Former", "Current"))
d$smoking_ever <- factor(d$smoking_ever, levels = c("Never", "Ever"))
d$dm           <- factor(d$dm)
d$htn          <- factor(d$htn)
d$statin       <- factor(d$statin)
d$circulation  <- factor(d$circulation)

## --- save (de-identified: no MRN/DOB/dates ever enter d) -----------------
dir.create("../outputs", showWarnings = FALSE)
saveRDS(d, "../outputs/analysis_data.rds")

## --- verification: reproduce manuscript Table 1 --------------------------
cat("\n================ Table 1 reproduction (target in brackets) ================\n")
cat(sprintf("N = %d  [171]\n", nrow(d)))
cat(sprintf("Age mean(sd): %.1f (%.1f)   [54 (15)]\n", mean(d$age, na.rm=TRUE), sd(d$age, na.rm=TRUE)))
cat("Sex:      "); print(table(d$sex, useNA="ifany")); cat("   [Female 113, Male 58]\n")
cat("Smoking:  "); print(table(d$smoking, useNA="ifany")); cat("   [Never 92, Former 22, Current 51, miss 6]\n")
cat("Diabetes: "); print(table(d$dm, useNA="ifany")); cat("   [No 147, Yes 19, miss 5]\n")
cat("HTN:      "); print(table(d$htn, useNA="ifany")); cat("   [No 49, Yes 84, miss 38]\n")
cat("Statin:   "); print(table(d$statin, useNA="ifany")); cat("   [No 125, Yes 40, miss 6]\n")
cat("HH miss=", sum(is.na(d$HH)), " mFS miss=", sum(is.na(d$mFS)),
    " WFNS miss=", sum(is.na(d$WFNS)), " GCS miss=", sum(is.na(d$GCS)), "  [38/38/35/39]\n")
cat(sprintf("MaxDiam mean(sd): %.1f (%.1f)  [7.1 (4.6)]\n", mean(d$maxdiam,na.rm=TRUE), sd(d$maxdiam,na.rm=TRUE)))
cat(sprintf("Aspect ratio mean(sd): %.2f (%.2f)  [2.02 (0.87)]\n", mean(d$aspect_ratio,na.rm=TRUE), sd(d$aspect_ratio,na.rm=TRUE)))
cat(sprintf("Size ratio mean(sd): %.2f (%.2f)  [2.05 (1.22)]\n", mean(d$size_ratio,na.rm=TRUE), sd(d$size_ratio,na.rm=TRUE)))
cat(sprintf("Volume mean(sd): %.0f (%.0f)  [257 (687)]\n", mean(d$volume,na.rm=TRUE), sd(d$volume,na.rm=TRUE)))
cat("Bleb present:", sum(d$bleb==1, na.rm=TRUE), " miss=", sum(is.na(d$bleb)), "\n")
cat("===========================================================================\n")
