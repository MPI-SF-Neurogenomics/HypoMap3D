# 11_Brainrender

This module is used for plotting [Brainrender](https://brainglobe.info/documentation/brainrender/index.html) figures.

Unfortunately, this part of the repo must be run on a local machine currently as we have not yet fixed technical issues that arose when trying to the scripts on our HPC (e.g., visual bugs in the plots etc.).

## Environment

This part is intended to be used without any of our containers.

The code was developed on Python 3.13 with the dependencies specified in requirements.txt.

## Running code

To run the code you need to

- Run 0_precompute_coord_space.py and 1_precompute_cells_per_region.ipynb on your HPC
- Clone the repo on your local machine 
- Set ../../config.env's DATADIR variable to a path on your machine where you want to save the objects donwloaded in the next step
- Download HypoMap3D's spatial cell-level and voxel-level datasets as well as the "coord_space.csv" and "coord_space_ext.csv" files created in the first step
- Place the spatial objects into DATADIR/9_objects/spatial/ and the .csv files into DATADIR/11_brainrender/
- __ensure that the current working directory is this directory__

### Running scripts specific to a figure

Scripts that are specific to a figure have their Figure in their name.

Such scripts can be run with

> python F1_2D_cell_count_heatmap.py 

### Running scripts unspecific to a figure

Scripts that are unspecific to a figure start with 2D_ or 3D_.

You can run these while passing a parameters file as an argument.

> python 2D_vox.py params/Ext_Fig9_arh_2D.json

The figures are using the following scripts:

Fig1:
- 3D_vox.py params/F1_agrp_vox_rm.json
- 3D_cells.py params/F1_agrp_cells.json

F3:
- 2D_vox.py params/F3_gsx1_2D.json

F5:
- 2D_vox.py params/F5_NN_2D.json

Ext Data Fig6:
- 2D_cells_and_vox_probs.py: all params/Ext_Fig6_ files with "cells_and_vox" in their name
- 2D_cells_two_genes.py: all params/Ext_Fig6_ files without "cells_and_vox" in their name

Ext Data Fig7-10:
- 2D_vox.py with their respective param files

Plots will be dumbed into DATADIR/11_brainrender/3D_plots or DATADIR/11_brainrender/2D_plots, respectively.
