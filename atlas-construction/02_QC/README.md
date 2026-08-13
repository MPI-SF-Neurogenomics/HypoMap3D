# 2_QC

In this module, Cellbender is run on all samples to remove empty droplets. 

Afterwards, standard quality control thresholds are used to remove low-quality cells.

Lastly, sample-wise h5ads are concatenated.

### Notes on usage

Please consider the docstrings in 1_CellBender.sh and 1a_submit_all_samples_cellbender.sh

In general, 1_CellBender.sh runs a process in your current terminal for a single sample while 1a will submit jobs for all samples to slurm.

Running either of the scripts suffices but using SLURM is heavily encouraged.

__Important: use the --cap-cells flag when running these scripts__