# CRC S09 Patient-Level Validation

This repository contains the frozen S09 reproducibility release for the study **Patient-level validation of colorectal cancer epithelial state programs across independent single-cell cohorts**.

The release documents a patient-level evaluation of colorectal cancer epithelial expression-state programs, a trajectory-execution audit, and separate descriptive progression analyses across public single-cell cohorts. Patients or donors are the inferential units; cells are used for coverage and quality-control accounting and are not treated as independent biological replicates.

## Scope

- GSE294300 is the independent paired validation cohort.
- GSE315534, GSE178318 and the Moorman 2024 cohort are separate descriptive progression cohorts.
- The analysis is frozen and reporting-only. No parameter tuning, pooled-cell inference or post hoc validation is included.
- No directional model, parameter tuning or post hoc validation is part of this release.

## Repository contents

- `code/`: deterministic figure and supplementary-table generation scripts.
- `config/frozen/`: frozen S08/S09 analysis contracts.
- `provenance/`: handoffs, freeze manifest, package status and audit records.
- `source_data/`: accession metadata, sample manifests and checksum records.
- `environment/`: execution-environment capture and package information.
- `tables/`: machine-readable supplementary tables and manifests.
- `reproduce.md`: reproducibility workflow and verification steps.
- `LICENSE`: release license.

## Data access

Large public count matrices and raw archives are not duplicated in this repository. Retrieve them from the original repositories using the accession identifiers and URLs in `source_data/`. Public-data access remains subject to the terms of the originating repositories.

## Reproducibility

The formal execution environment was WSL Ubuntu-24.04 using the project `cast-py` and `cast-r` environments on CPU. See `reproduce.md` for the frozen workflow, input boundaries and verification steps. The rendering scripts require the public inputs and the frozen analysis handoff files described in the provenance records. SHA-256 records are provided where available.

## Release status

This is a frozen S09 reporting release. It is intended to support manuscript review and independent inspection of the reported figures, tables and analysis contracts.
