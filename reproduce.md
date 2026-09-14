# Reproducibility Workflow

## 1. Environment

Use Ubuntu-24.04 under WSL2 with the project environments specified in `environment/`. The formal analysis used CPU execution. Do not substitute another project environment without recording the change.

## 2. Retrieve public inputs

Read the accession identifiers, sample boundaries, source URLs and checksum records in `source_data/`. Download large matrices only from their original repositories. Do not infer patient boundaries from cell barcodes without consulting the corresponding manifest.

## 3. Inspect frozen contracts

Review the applicable files in `config/frozen/` before running any script. The contracts define cohort roles, input boundaries, analysis units, minimum criteria and interpretation restrictions.

## 4. Regenerate reporting outputs

After the public inputs and project environments are available, run the scripts in `code/` with the documented project environment. The scripts generate reporting figures and supplementary tables; they do not retune the frozen analysis.

## 5. Verify the release

Compare generated files with the provenance records and verify available SHA-256 values. Preserve failed or abstained trajectory executions in the audit rather than removing them.

## Interpretation boundary

Expression-state summaries, trajectory execution outputs and progression-associated estimates are separate evidence layers. A finite pseudotime output is an execution result and does not by itself establish lineage, conversion probability or biological direction.
