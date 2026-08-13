# 08 Annotation

In this module, we annotate our clusters with cell type labels.

For neuronal cells, we provide two different versions: ABC-like and descriptive short-form annotations. For details, please see our manuscript.

Afterwards, we finish by calculating region-specific marker genes.

### Notes on script execution

3_submit_region_specific_markers.sh will submit SLURM jobs for the calculation region-specific markers for all levels and the regions specified within that script, which are all regions of interest for us (hypothalamic regions or regions that were used for short-form annotation of L3 cell types).

Jobs will use the 3_region_specific_markers.py script, passing levels and regions as arguments. Hence, using the SLURM wrapper suffices but you may also use the .py script without the SLURM wrapper if you want to calculate markers for a specific region or level.