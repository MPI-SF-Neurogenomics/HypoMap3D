# 04_HY

In this module, we aim to remove cells that are likely to come from non-hypothalamic regions.

For that matter, we first train a new scVI model, after which Allen Brain Cell cell types are propagated to other cells.

Using the regional annotation within the ABC cell type names, we remove cells likely being non-hypothalamic.


### Notes on usage

For the scVI part, the script with "submit" in its name will submit a job to SLURM while the other one will run scVI in the terminal. 

Running either of the scripts suffices but using SLURM is heavily encouraged.