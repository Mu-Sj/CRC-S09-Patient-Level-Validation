# Submission Package Consistency Audit

Audit date: 2026-09-13

Scope: read-only cross-check of the current `submission_package/` against its manuscript, legends, tables, frozen handoffs, active-file index and package checksums. No analysis, parameter, cohort role or figure was changed during this audit.

## Findings

### Resolved: trajectory replicate-count interpretation

The frozen GSE294300 final handoff reports 11 total execution records (replicate 0 plus 10 patient-bootstrap replicates). `tables/SupplementaryTable4.tsv` reports `replicates=10.0` for PAGA, Palantir, Slingshot and Monocle3 because that field refers to the 10 bootstrap attempts. The manuscript, unified figure legend and supplementary methods now state both quantities explicitly. The 10-of-10 readiness counts refer to bootstrap attempts; the 11-record denominator is used only for the complete execution audit. This interpretation agrees with `run_metadata.txt` (`replicates=0_plus_10`) and the 11 rows in the formal trajectory records.

### Medium priority: stale historical audit text

`provenance/S12_SUBMISSION_GAP_AUDIT_v2.md` still states that the manuscript, figures, supplementary materials, release package and submission statements are missing. Those materials now exist in the package. The document is an historical audit; `S12_CURRENT_PACKAGE_STATUS.md` prevents the old `MISSING` entries from being read as current.

`provenance/S12_FIGURE_STORYBOARD.md` still describes a four-figure partial set and says that the files are not final submission figures. The active index now specifies five main figures and six supplementary figures. The storyboard is retained as design history; `S12_CURRENT_PACKAGE_STATUS.md` identifies the active set.

### Medium priority: active-file index scope

`CURRENT_SUBMISSION_FILES.tsv` now indexes the five active main PDFs, six supplementary PDFs, the table manifest, the three reproducibility scripts, the manuscript files, supplementary text and submission statements. It is the complete current journal-facing file manifest; raw source matrices remain outside the package.

### Passed checks

All active paths listed in `CURRENT_SUBMISSION_FILES.tsv` exist. The package SHA-256 manifest passes. The package contains no file larger than 100 MB and no raw count matrix or raw archive. The manuscript, unified figure legends, supplementary methods and statements use the patient/donor inferential-unit rule and do not introduce a pooled-cell result or a new analysis. The primary endpoint values agree across the manuscript, Figure 3 legend, supplementary reporting notes and the frozen GSE294300 handoff: absorptive mean difference -0.20474 (95% CI -0.26994 to -0.14409) and secretory mean difference +0.05494 (95% CI -0.02083 to +0.15804).

## Required resolution before release

1. Confirm whether the authoritative formal trajectory denominator is 11, as stated in the final handoff, or 10, as stated in `SupplementaryTable4.tsv` and the current manuscript readiness paragraph.
2. Reconcile `SupplementaryTable4.tsv`, manuscript text, Figure 4 legend and Supplementary Figure 4 legend to the confirmed denominator. Recompute any displayed percentages only from the confirmed frozen record; do not retune or rerun the analysis.
3. Retain `S12_CURRENT_PACKAGE_STATUS.md` alongside the historical gap audit and storyboard so that their old `MISSING` and four-figure statements cannot be read as current.
4. Keep `CURRENT_SUBMISSION_FILES.tsv` as the complete journal-facing file manifest.

## Audit conclusion

The package is structurally complete and checksum-valid. The trajectory count has been reconciled as 10 bootstrap replicates plus 1 initial record, a current status note supersedes historical missing-item language, and the active index now includes all journal-facing text and statement files. No new analysis is needed; only author-specific fields and target-journal formatting remain.
