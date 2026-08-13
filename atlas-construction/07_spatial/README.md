# 07 Spatial

In this module, we leverage the datasets listed below to generate spatial predictions for our voxel space.

In short, we source and harmonize the datasets after which we transfer labels from our single cell atlas.

Then, cells are assigned to voxels spanning the hypothalamus boundary according to ABC's Common Coordinate Framework v3, and the cell's cell type membership probabilities aggregated per voxel. Probabilities are then smoothed in 3D, corrected for cell densities and normalized, resulting in voxel-wise relative masses for all cell types.

Lastly, we analyze for region enrichment intended for cell type annotation and derive metrics to judge the fidelity of our predictions, culminating in our spatial warnings/flags.

## Datasets

### SP2-01 (Yao | MERFISH | ABC)

> Publication: [Paper Link](https://www.nature.com/articles/s41586-023-06812-z)  
> Data Access: [Data Link](https://alleninstitute.github.io/abc_atlas_access/notebooks/merfish_tutorial_part_2a.html)

### SP2-02 (Zhang | MERFISH | ABC)

> Publication: [Paper Link](https://www.nature.com/articles/s41586-023-06808-9)  
> Data Access: [Data Link](https://alleninstitute.github.io/abc_atlas_access/notebooks/zhuang_merfish_tutorial.html)

### SP2-03 (Bauder | OpenST | MPI)

> No publication available. Dataset generated in-house.  
> Data access will be made available on GEO.

### SP2-04 (Langlieb | SlideSeq | BI)

> Publication: [Paper Link](https://www.nature.com/articles/s41586-023-06818-7)  
> Data Access: [Data Link](https://drive.google.com/drive/folders/1ajep-hCLW1BJ-ljS_3ApjgpVV4uELpFg)

