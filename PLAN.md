# Resubmission Execution Plan

Approved sequence: **(1) Tables + Figures → (2) Manuscript → (3) Response letter.**
Autonomous execution once approved. Every artifact is regenerated from the raw
data by the R pipeline (reproducible). Relevant skills are applied per stage.

## Decision log
- **Sparse PCA / PLS / clustering — REJECTED** (critical-thinking + senior-DS review).
  Sparse PCA is still unsupervised (won't beat PC1); PLS is circular without nested
  CV (overfits to the outcome); clustering on 7 size-dominated indices just re-derives
  "big vs small." All three add multiplicity + a method-shopping appearance without
  advancing the thesis. Instead: keep the parsimonious pre-specified size-ratio/aspect-
  ratio ordinal model, and *explain* PC1's behavior (variance ≠ association) in the
  Discussion; note PC1 IS ordinally associated with WFNS/GCS in the response letter.

## Stage 0 — Done (foundation)
- ✅ Docs converted (pandoc); reviewer comments analyzed.
- ✅ `01_clean.R` reproduces manuscript Table 1 exactly (validated recoding).
- ✅ `02_analysis.R` corrected + extended models (ordinal + MICE + PCA + threshold + index).
- ✅ `03_rwe_checks.R` MICE bias audit (MI-then-delete, MNAR delta, E-values) — no bias.
- ✅ `docs/METHODS_AND_EQUATIONS.md`, `docs/CODE_WALKTHROUGH.md`.
- ✅ Figures 1–3 drafted and rendering at 300 dpi.

## Stage 1 — Tables + Figures  (skills: scientific-visualization, statistical-analysis, visualization-best-practices, accessibility)
1. **`05_tables.R`** — publication tables, all regenerated from data:
   - **Table 1** — cohort descriptives (already validated).
   - **Table 2 (corrected)** — adjusted proportional-odds ordinal ORs (per SD),
     complete-case + MICE, with FDR. Replaces the mislabeled "ridge (LASSO)" table.
   - **Table 3 (corrected)** — indices surviving FDR under MICE (the honest headline).
   - **Supp tables** — full grid, PC1/PC2, logistic dichotomized, E-values, sensitivity.
   - Output: CSV + Word-ready `.docx` via a small render step.
2. **Figures** — finalize the three built (fix em-dash glyphs → clean TIFF),
   colorblind-safe, 300 dpi TIFF + vector PDF. Confirm each visually.
3. **Method figure (Fig 4)** — analytic-workflow schematic (cohort/exclusions →
   7 morphometrics → ordinal + MICE → sensitivity), for manuscript + repo.
   (skill: scientific-schematics / diagrams).
4. Write `figures/FIGURE_LEGENDS.md`.

## Stage 2 — Repo finalize  (skills: python-project-structure/repro, audit-reproducibility)
- Wire `04_figures.R`/`05_tables.R` into a `run_all.R`; add `renv.lock`;
  update README with the figure/table map. Keep PHI-safe (.gitignore).

## Stage 3 — Manuscript rewrite  (skills: scientific-writing, the-humanizer [formal journal tone], storytelling, scientific-critical-thinking)
- **Abstract** — reframe to the corrected positive finding (ordinal + MICE);
  fix numbers; drop the internal contradiction.
- **Methods** — ordinal regression + PO check; MICE (primary) + complete-case
  (sensitivity); standardized predictors; ridge for ranking only; E-values; STROBE.
- **Results** — new ordinal/MICE ORs; PC1 size-axis vs PC2 shape-axis; threshold;
  honest negatives (bleb null, AUC ~0.56); reference new Tables/Figures.
- **Discussion** — size-ratio-dominant mechanism; reconcile the composite;
  unmeasured-confounding via E-values; corrected limitations.
- **Global fixes** — "ASPECTS ratio" → **aspect ratio** throughout; correct Table 2
  method label; remove invalid ridge p-values.
- Deliver as tracked-changes-ready `.docx`.

## Stage 4 — Response-to-reviewers letter  (skills: respond-to-referees, peer-review, scientific-writing)
- Point-by-point to R1 (contradiction, missingness, table discrepancy, ridge/LASSO)
  and R2 (ordinal regression, assumption checks), each citing the new numbers and
  pointing to the revised table/figure. Professional, concise, non-defensive tone.

## Guardrails during autonomous run
- Raw PHI xlsx never committed or moved; all outputs de-identified.
- Every reported number traces to a script; no hand-typed statistics.
- I will NOT submit anything, push to GitHub, or email — those remain your action.
- I pause and surface if any result materially contradicts the current narrative.
