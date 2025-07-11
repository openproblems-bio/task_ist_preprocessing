

Here we generate a small test dataset, used for `viash test`. Note that the file structure here is a bit simplified compared to `scripts/create_resources` as we only have one dataset.

Download and process the single cell data:
`bash 2023_yao_mouse_brain_scrnaseq_10xv2.sh`

Download and process the spatial data:
`bash 2023_10x_mouse_brain_xenium_rep1.sh`

Combine the two datasets and run the ist preprocessing pipeline once with generic methods to create example outputs after each step: `test_pipeline.sh`
