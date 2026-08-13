from pathlib import Path
import os 
import pandas as pd 
import brainglobe_heatmap as bgh
import brainrender
import numpy as np
import matplotlib.pyplot as plt 
from brainrender.actors import Points,PointsDensity
from dotenv import dotenv_values

brainrender.settings.SHADER_STYLE = "cartoon" 

#################################
####### only showing the HY structures here! this can be misleading
#################################

ANNOTATE_REGIONS = False
SLICE_THICKNESS = 250
PLAY = int(SLICE_THICKNESS / 2)
distances = list(range(4750, 8500, SLICE_THICKNESS))
print(distances)
# Define desired x/y limits
XLIM = (-100, 2600)  # only right
YLIM = (3200, -500)

params = {
# SCH surpasses vmax at 6000 with 1.9893048128342246
# ARH surpasses vmax at 6500 with 1.8106701366297984
# MMme surpasses vmax at 8250 with 2.0252100840336134
    "st": {
        "cell_counts_column": "ST_cell_counts",
        "cmap": "viridis",
        "vmax": 1.8, # when aggregated per brain region, SCH is 1.55, second one is 1.3..
        "suptitle": "ST cells along the AP axis",
        "suptitle_fontsize": 24

    },

# SCH surpasses vmax at 6000 with 6.483515559632937
    "sc": {
        "cell_counts_column": "estimated_sc_cell_counts",
        "cmap": "plasma",
        "vmax": 5, # when aggregated per brain region, SCH is at 5.6, ARH at 4.2; 
        "suptitle": "Estimated sc/sn cells along the AP axis",
        "suptitle_fontsize": 22
    }
}

CONFIG = dotenv_values("../../config.env") # assuming you are in the script on local machine without uv container
DATADIR = Path(CONFIG['DATADIR'])

SPATIALDIR = DATADIR / "9_objects" / "spatial"
BRAINRENDER_DIR = DATADIR / "11_brainrender" 
OUTDIR = BRAINRENDER_DIR / "2D_plots"
os.makedirs(OUTDIR, exist_ok=True)

df = pd.read_csv(BRAINRENDER_DIR / "coord_space_ext.csv")

def get_values_for_slice(spatial_slice, column_name):
    num_vox_per_nuc = spatial_slice.groupby("allen_acronym")[column_name].count()
    num_vox_per_nuc.name = "num_vox_per_nuc"
    num_cells_per_nuc = spatial_slice.groupby("allen_acronym")[column_name].sum()
    num_cells_per_nuc.name = "num_cells_per_nuc"
    scatter_df = pd.concat([num_vox_per_nuc, num_cells_per_nuc],axis=1).reset_index()
    scatter_df.allen_acronym = scatter_df.allen_acronym.astype(str)
    scatter_df["cells_per_vox"] = scatter_df.num_cells_per_nuc / scatter_df.num_vox_per_nuc
    return scatter_df.set_index("allen_acronym").cells_per_vox.to_dict()

for k, param_dict in params.items():

    print(f"Processing {k}")
        
    scenes = []
    vmax = param_dict["vmax"]
    cmap = param_dict["cmap"]
    column_name = param_dict["cell_counts_column"]
    suptitle = param_dict["suptitle"]
    suptitle_fontsize = param_dict["suptitle_fontsize"]

    for distance in distances:

        #get a slice
        spatial_slice = df.loc[(df.x_ccf > distance - PLAY) & (df.x_ccf < distance + PLAY) & (
        df.is_HY_structure),: ]
        
        values = get_values_for_slice(spatial_slice, column_name=column_name)

        # Step 1: Sort by value (descending)
        top_items = sorted(values.items(), key=lambda x: x[1], reverse=True)[:5]
        # Step 2: Print the top 2 key–value pairs
        for key, value in top_items:
            if value > vmax:
                print(f"{key} surpasses vmax at {distance} with {value}")

        scene = bgh.Heatmap(
            values,
            position=distance,
            orientation="frontal",
            thickness=SLICE_THICKNESS,
            format="2D",
            cmap=cmap,
            vmin=0,
            vmax=vmax,
            annotate_regions=ANNOTATE_REGIONS,
            hemisphere="left",
            root=False,
        )

        scenes.append(scene)

    # Create a figure with 6 subplots and plot the scenes
    num_cols = int(np.ceil(len(distances)/3))
    num_rows = 3
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(num_cols * 1.4, num_rows * 2))

    for scene, ax in zip(scenes, axs.flatten(), strict=False):
        scene.plot_subplot(fig=fig, ax=ax, show_cbar=False, hide_axes=True)
        
        # Set axis limits AFTER plotting
        ax.set_xlim(XLIM)
        ax.set_ylim(YLIM)
        
        ax.set_title(ax.get_title(), fontsize=13, 
        ) 

    fig.suptitle(suptitle, fontsize=suptitle_fontsize, weight='bold'
    )

    plt.tight_layout()

    fig_save_path = OUTDIR / f"{k}_cell_density_2D_across_brain_heatmap.pdf"
    fig.savefig(fig_save_path, dpi=300)

    plt.show()

