

Here we generate a small test dataset, used for `viash test`. Note that the file structure here is a bit simplified compared to `scripts/create_resources` as we only have one dataset.

Download and process the single cell data:
`bash 2023_yao_mouse_brain_scrnaseq_10xv2.sh`

Download and process the spatial data:
`bash 2023_10x_mouse_brain_xenium_rep1.sh`

Combine the two datasets and run the ist preprocessing pipeline once with generic methods to create example outputs after each step: `test_pipeline.sh`

Two additional spatial test samples (cropped to a small region and published under `resources_test/common/`):
- 10x Atera WTA FFPE human breast cancer: `bash 2026_10x_human_breast_cancer_atera.sh`
- Bruker CosMx human lung cancer (Lung9 Rep1, smallest CosMx archive): `bash 2026_bruker_human_lung_cancer_cosmx.sh`

Each of the two also has a `*_nebius.sh` variant that runs on Nebius via `tw launch` and publishes to `/scratch/task_ist_preprocessing/resources_test/common/` instead of running locally (recommended for cosmx, which is too heavy for a laptop). After the run, sync the result from scratch to S3 (see the commented `aws s3 sync` at the bottom of each script).

Matching single-cell references for the two spatial samples (Nebius `tw launch` → scratch; each script runs the SC workflow, then chains the `subsample` module down to ~400 cells since these workflows don't subsample internally):
- Breast cancer reference for atera: `bash 2021_wu_human_breast_cancer_scrnaseq_nebius.sh`
- Lung cancer (NSCLC) reference for cosmx: `bash 2024_zuani_human_nsclc_scrnaseq_nebius.sh`

Pairings for a combined dataset (à la `mouse_brain_combined`): atera ↔ Wu breast, cosmx lung ↔ Zuani NSCLC.
