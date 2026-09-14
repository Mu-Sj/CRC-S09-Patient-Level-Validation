# S09 Safe-Path Handoff

Date: 2026-09-13  
Status: `CONDITIONAL_DESCRIPTIVE_ONLY`

## Scope

The formal S09 primary analysis remains the frozen patient/donor-level state-program endpoint. For each eligible patient or donor, the mean module score is computed within each arm using CP10K normalization, `log1p`, and the frozen S05 genes, followed by a paired arm difference and complete-patient bootstrap percentile interval.

The progression cohorts are retained as separate descriptive analyses. Their estimates are not pooled because treatment, sampling, malignant-cell evidence, and processing differ across cohorts.

## Technical sensitivity method

Monocle3 is the single predeclared technical-sensitivity method. This is supported by 10/10 formal S08 patient-bootstrap replicates and 11/11 S09 GSE294300 replicates passing the declared finite-pseudotime checks. The method is used only to report execution stability, root dependence, and failure/abstention records.

Monocle3 is not used to select a biological direction, retune the S05 modules, change thresholds, replace the patient-level primary endpoint, or claim lineage conversion.

## Results retained for reporting

| Cohort | Edge | Module | n | Mean delta | Bootstrap 95% CI | Status |
|---|---|---|---:|---:|---:|---|
| GSE315534 | primary -> liver metastasis | absorptive | 3 | -0.1134 | -0.1595 to -0.0891 | descriptive |
| GSE315534 | primary -> liver metastasis | secretory | 3 | -0.1946 | -0.2184 to -0.1772 | descriptive |
| GSE178318 | primary -> liver metastasis | absorptive | 6 | -0.0101 | -0.0394 to +0.0187 | descriptive |
| GSE178318 | primary -> liver metastasis | secretory | 6 | -0.0173 | -0.0845 to +0.0417 | descriptive |
| Moorman 2024 | primary -> metastasis | absorptive | 6 | +0.0148 | -0.0284 to +0.0543 | descriptive |
| Moorman 2024 | primary -> metastasis | secretory | 6 | -0.0297 | -0.0846 to +0.0179 | descriptive |

## Manuscript positioning

The defensible claim is that frozen epithelial state programs can be measured reproducibly at the patient/donor level, while progression deltas show cross-cohort non-replication and substantial uncertainty. Monocle3 contributes a method-execution sensitivity analysis, not a proof of a metastasis lineage or conversion probability.

The manuscript should explicitly state the small number of independent patients, the restricted/heterogeneous malignancy evidence, the GSE178318 COL17 low-count metastasis sensitivity case, mixed treatment strata, and the fact that S09 is not a formal progression holdout.

## Prohibited changes

No pooled-cell p-values, post hoc patient removal, S05 module changes, threshold changes, validation-result-driven root selection, software ranking, composite score, velocity/CellRank activation, or S10 activation is authorized.

## Reproducibility anchors

- Primary S09 handoff: `analysis/09_external_validation/S09_GSE294300_FINAL_HANDOFF.md`
- Cross-dataset summary: `analysis/09_external_validation/runs/20260913_S09_direction_audit_v1/cross_dataset_primary_to_metastasis_summary.tsv`
- Direction audit: `analysis/09_external_validation/runs/20260913_S09_direction_audit_v1/S09_DIRECTION_AUDIT.md`
- Monocle3 decision: `docs/decisions/20260913_S09_monocle3_method_selection.md`
