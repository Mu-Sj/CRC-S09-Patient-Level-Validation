# S09 GSE294300 Final Handoff and Readiness

Date: 2026-09-12  
Dataset: GSE294300  
Status: `CONDITIONAL_DESCRIPTIVE_ONLY`

## Decision

The frozen S09 primary route is complete for a patient-level paired state
concordance endpoint and a secondary trajectory execution audit. The cohort is
ready for descriptive handoff only. It is not ready for a real-CRC software
ranking, a composite score, a lineage claim, or a conversion-probability claim.

## Frozen scope and independence

- Primary denominator: 16 patients with paired adjacent-normal and primary-tumor
  samples.
- Sensitivity layer: patients 1 and 17; retained outside the primary denominator
  because the frozen sample-level candidate floor was not met.
- Contract: `config/frozen/S09_GSE294300_validation_contract_v1.1_20260912.yaml`.
- Input: `data/processed/GSE294300_lenient_candidates_v1`; raw counts remain
  preserved separately.
- Patient is the inference unit; cells are observations within patients.
- GSE294300 and the GSE236581 discovery cohort have distinct GEO study/sample
  accessions and no literal patient-key overlap in audited metadata. Public
  metadata cannot prove that no participant contributed to both studies.

## Primary endpoint: paired frozen-state concordance

The tumor-minus-normal patient-level module-score difference used 1,000 complete-
patient bootstrap resamples.

| Frozen module | Mean paired difference | Patient-bootstrap 95% CI |
|---|---:|---:|
| `SECRETORY_MUC2_TFF3` | +0.05494 | -0.02083 to +0.15804 |
| `ABSORPTIVE_BEST4_CA4` | -0.20474 | -0.26994 to -0.14409 |

The absorptive program decrease is precise under this contract; the secretory
increase is positive but imprecise because its interval includes zero. The result
supports reproducibility of state programs in this cohort, not a biological
lineage transition.

## Secondary endpoint: trajectory execution and patient bootstrap

All trajectory runs used the frozen Normal-root rule, complete-patient bootstrap,
per-sample cap of 2,000 cells, 2,000 HVGs, 15 PCA dimensions and no dense count
conversion.

| Runner/method layer | Replicates | PASS | ABSTENTION | INVALID_RESULT | HARD_FAILURE | Notes |
|---|---:|---:|---:|---:|---:|---|
| R Slingshot | 11 | 11 | 0 | 0 | 0 | 3--6 lineages; union of finite pseudotimes covered all cells |
| R Monocle3 | 11 | 11 | 0 | 0 | 0 | finite pseudotime for every cell; graph metrics retained |
| Python PAGA + Palantir runner | 11 | 9 | 0 | 0 | 2 | replicate 0 and 8 failed at Palantir/ARPACK eigensolver convergence |

The R audit is `analysis/09_external_validation/runs/20260912_S09_gse294300_r_trajectory_formal_v9/S09_R_TRAJECTORY_AUDIT.md`.
The Python records are retained in
`analysis/09_external_validation/runs/20260912_S09_gse294300_python_trajectory_formal_v1/trajectory_pilot.tsv`
and `summary.json`; the two ARPACK failures are part of the result, not removed
or repaired after inspection.

These outputs establish execution and patient-bootstrap stability for the tested
layers. They do not establish direction, lineage truth, conversion probability,
or a software winner. No real-CRC composite score or ranking is generated.

## Quality controls and limitations

- The analyzed matrix is a provisional barcode-rank-derived candidate input; it
  is not asserted to reproduce an author-provided filtered matrix.
- Ambient-RNA correction and final doublet adjudication remain missing/review
  evidence in this route.
- CNV is auxiliary malignancy/genomic-consistency evidence and is not malignant
  cell truth, a state label, or a direction label.
- Python trajectory execution had an 18.2% hard-failure rate (2/11) from ARPACK
  non-convergence. The frozen failure contract requires preservation of these
  records; no validation-result-driven parameter change was made.
- R logs contain non-fatal Slingshot Quick-TRANSfer warnings; all 11 R replicates
  still passed the declared finite-coverage checks.
- Public metadata support cohort-level independence but cannot exclude real-world
  participant overlap across studies.

## Reproducibility record

Environment: WSL `Ubuntu-24.04`; Python `/home/msj_proj/miniconda3/envs/cast-py`;
R `/home/msj_proj/miniconda3/envs/cast-r`; CPU execution; no GPU; formal R
bootstrap used the complete-patient unit and 11 replicates (replicate 0 + 10).

Formal R command:

```bash
/home/msj_proj/miniconda3/bin/conda run --no-capture-output \
  -p /home/msj_proj/miniconda3/envs/cast-r Rscript \
  /mnt/e/bis_project/CAST/analysis/09_external_validation/scripts/gse294300_r_runner.R
```

R run metadata: `analysis/09_external_validation/runs/20260912_S09_gse294300_r_trajectory_formal_v9/run_metadata.txt`  
R log: `logs/20260912_S09_gse294300_r_trajectory_formal_v9.log`  
R exit code: `0`  
State-endpoint log: `logs/20260912_S09_gse294300_state_concordance.log`

Key package versions: Python `anndata 0.12.16`, `scanpy 1.12.1`, `scvelo 0.3.4`,
`palantir 1.4.5`, `numpy 2.3.5`, `scipy 1.16.3`; R
`SingleCellExperiment 1.32.0`, `scran 1.38.0`, `scater 1.38.0`,
`slingshot 2.18.0`, `monocle3 1.4.27`, `Matrix 1.7.5`.

SHA-256 (key machine-readable outputs):

```text
987F2BD279800A1150519CA8799DFDE11D3CF92771BBF527B8947A5FD819B860  state_concordance/summary.json
89AF476FD7F8E66D085EF439BAA8F5000B61EA4EEF733A43DE025038E6C90E17  state_concordance/patient_paired_state_concordance.tsv
64B500DAA91A67311714AEE686391F6FDA842E8FDC505597CCDFE09045549B0A  r_trajectory_formal_v9/r_trajectory.tsv
118698BA63A53B3F7F1E2FC912B1A0A8C0CE90B7DB4283F9B0EA8952A7A2CC75  python_trajectory_formal_v1/trajectory_pilot.tsv
```

## Next gate

S09 can be marked complete for this descriptive handoff. Keep the project state
`CONDITIONAL_DESCRIPTIVE_ONLY`. Any progression analysis using GSE178318 needs
its own epithelial/CNV evidence audit and progression-specific frozen contract;
it is not part of this handoff. No further S09 tuning or ranking should start
without a newly approved scope and an independent contract.
