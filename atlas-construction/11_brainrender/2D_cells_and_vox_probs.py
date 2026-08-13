
from pathlib import Path
import os 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable, get_cmap
from matplotlib.patches import Polygon
import anndata as ad
import scanpy as sc
import argparse
import zarr
import json
import brainglobe_heatmap as bgh
from dotenv import dotenv_values
from matplotlib.ticker import ScalarFormatter

CONFIG = dotenv_values("../../config.env") # assuming you are in the script on local machine without uv container
DATADIR = Path(CONFIG['DATADIR'])

SPATIALDIR = DATADIR / "9_objects" / "spatial"
BRAINRENDER_DIR = DATADIR / "11_brainrender" 
OUTDIR = BRAINRENDER_DIR / "2D_plots"
os.makedirs(OUTDIR, exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument(
    "params",
    help="Path to parameters JSON file."
)

args = parser.parse_args()

with open(args.params) as f:
    params = json.load(f)
    
SLICE_THICKNESS = params["SLICE_THICKNESS"] #step size from start to stop range
START_RANGE = params["START_RANGE"] # 4750 is default, 6200 for periARH, 6500 for tany
STOP_RANGE = params["STOP_RANGE"] # 8500 is default, 7800 for periARH, 7800 for tany
PLAY_FOR_CELLS = params["PLAY_FOR_CELLS"] # margin for cells to be plotted on the slice
ANNOTATE_REGIONS = params["ANNOTATE_REGIONS"] # False is default
SCATTERSIZE_CELL = params["SCATTERSIZE_CELL"] # 0.2 is default
SCATTERSIZE_VOX = params["SCATTERSIZE_VOX"] 
LINEWIDTH = params["LINEWIDTH"] #default 0.6
XLIM = params["XLIM"] # default [-1900, 1900]
YLIM = params["YLIM"] # default [3200, 800]
PROB_CMAP = params["PROB_CMAP"]
MIN_CELL_P_TO_PLOT = params["MIN_CELL_P_TO_PLOT"] # exclusive: cells must have p greater than this
CLUSTER_TO_PLOT = params["CLUSTER_TO_PLOT"]
MODALITY_VOX = params["MODALITY_VOX"]

# derive settings
ID = CLUSTER_TO_PLOT.split(" ")[0] # failsafe for if descr is passed
level_key = f"L{ID[2]}"
distances = list(range(START_RANGE, STOP_RANGE, SLICE_THICKNESS))
n_slices = len(distances)

# --- Load Data ---#
voxel_store = zarr.open(SPATIALDIR / 'hypomap3d_spatial_voxel_level.zarr')

adata = ad.read_h5ad(SPATIALDIR / "hypomap3d_spatial_cell_level.h5ad", backed="r")
celltype_idx = list(adata.uns[f"{level_key}_ID"]).index(ID)
cluster_descr = adata.uns[f"{level_key}_descr"][celltype_idx]

df = pd.read_csv(BRAINRENDER_DIR / "coord_space.csv")
df.allen_acronym = df.allen_acronym.fillna("HY") # this is a dummy such that regions not in the voxel_store will be included in the heatmap serving as a background
vx_idx = voxel_store[level_key]["ids"][:].tolist().index(ID)
probs_vox = voxel_store[level_key][MODALITY_VOX][vx_idx, :, :, :].reshape(-1)
probs_vox = np.nan_to_num(probs_vox, nan=0.0)
df[ID] = probs_vox
df = df.loc[df.z_ccf < 5700,:] # remove protruding voxels from right hemisphere

# extract data into .obs for convenience
point_array = adata.obsm["coords_ccf_mm_wrap"]
point_array = np.array(point_array, dtype = "float32")
point_array = point_array * 1000 
adata.obs["x_ccf"] = point_array[:, 0]
adata.obs["y_ccf"] = point_array[:, 1]
adata.obs["z_ccf"] = point_array[:, 2]
adata.obs[f"{ID}_prob"] = adata.obsm[f"P_{level_key}"][:, celltype_idx].toarray()

# filter data to only the slices to be plotted:
adata_mask = np.asarray([False] * len(adata))
df_mask = np.asarray([False] * len(df))
for x in distances:
    closest_x = df.x_ccf.unique()[np.abs(df.x_ccf.unique() - x).argmin()]
    min_x = closest_x - PLAY_FOR_CELLS
    max_x = closest_x + PLAY_FOR_CELLS
    adata_mask |= ((adata.obs.x_ccf > min_x) & (adata.obs.x_ccf < max_x))
    df_mask |= ((df.x_ccf > min_x) & (df.x_ccf < max_x))

adata = adata[adata_mask, "Tbx3"].to_memory() #just a dummy gene 
df = df.loc[df_mask, :]

# --- Get colormaps while getting vmax from all cells ---
def create_colormap(data, cmap_name):
    vmax = data.max()
    norm = Normalize(vmin=0, vmax=vmax)
    cmap = get_cmap(cmap_name)
    sm = ScalarMappable(norm=norm, cmap=cmap)
    return cmap, norm, sm

prob_colormaps = {}

prob_colormaps["cells"] = create_colormap(adata.obs[f"{ID}_prob"].values, PROB_CMAP)
prob_colormaps["vox"] = create_colormap(df[ID], PROB_CMAP)

# --- Precompute Heatmap Scenes ---
scenes_by_x = {}
for x in distances:
    closest_x = df.x_ccf.unique()[np.abs(df.x_ccf.unique() - x).argmin()]
    slice_data = df[df.x_ccf == closest_x]

    # Initialize with dummy values
    values = {k: 0 for k in slice_data.allen_acronym.unique()}
    
    scene = bgh.Heatmap(
        values,
        position=closest_x,
        orientation="frontal",
        thickness=1,
        format="2D",
        cmap="Greys",
        vmin=0,
        vmax=1,
        annotate_regions=ANNOTATE_REGIONS,
        hemisphere="both",
        root=False,
    )

    scenes_by_x[closest_x] = scene

# --- Plotting ---
fig = plt.figure(figsize=(n_slices * 4, 2.25))
height_ratios = [1.0, 0.2]
gs = GridSpec(2, n_slices + 3, figure=fig, 
    width_ratios=[1]*n_slices + [0.05, 0.02, 0.05], wspace=0.3,
    hspace=0.01, height_ratios=height_ratios)

# extract colormaps
cmap_prob_cell, norm_prob_cell, sm_prob_cell = prob_colormaps["cells"]
cmap_prob_vox, norm_prob_vox, sm_prob_vox = prob_colormaps["vox"]

for col, x in enumerate(distances):
    closest_x = df.x_ccf.unique()[np.abs(df.x_ccf.unique() - x).argmin()]
    scene = scenes_by_x[closest_x]
    ax = fig.add_subplot(gs[0, col])
    scene.plot_subplot(fig=fig, ax=ax, show_cbar=False, hide_axes=True)

    ax.set_xlim(XLIM)
    ax.set_ylim(YLIM)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')

    # --- RIGHT Hemisphere: (left on plot) cell probs ---
    prob_slc = adata[(adata.obs.x_ccf > closest_x - PLAY_FOR_CELLS) & (adata.obs.x_ccf < closest_x + PLAY_FOR_CELLS)]
    if prob_slc.shape[0] > 0:
        # ALIGN WITH HEATMAP
        _, y_bias, z_bias = scene.slicer.plane0.center_of_mass()
        coords = prob_slc.obs[["z_ccf", "y_ccf"]].values.copy()
        coords[:, 0] += z_bias
        coords[:, 1] -= y_bias

        # restrict to MIN_CELL_P_TO_PLOT
        prob_values = prob_slc.obs[f"{ID}_prob"].values
        min_p_mask = prob_values > MIN_CELL_P_TO_PLOT
        prob_values = prob_values[min_p_mask]
        coords_p = coords[min_p_mask]

        # sort so highest probability is plotted last
        order = np.argsort(prob_values)
        prob_values = prob_values[order]
        coords_p = coords_p[order]

        colors = cmap_prob_cell(norm_prob_cell(prob_values))
        ax.scatter(coords_p[:, 0], coords_p[:, 1], c=colors, s=SCATTERSIZE_CELL, alpha=1)

        # --- LEFT Hemisphere: (right side of plot) Cluster MODALITY (voxels) ---
        slc = df[(df.x_ccf == closest_x)]
        if not slc.empty:
            coords = slc[["z_ccf", "y_ccf"]].values.copy()
            colors = cmap_prob_vox(norm_prob_vox(slc[ID].values))
            # ALIGN WITH HEATMAP
            coords[:, 0] += z_bias
            coords[:, 1] -= y_bias
            coords[:, 0] = -coords[:, 0]  # LEFT hemisphere flip

            #print(f"[{cluster}] LEFT @ x={closest_x}: {len(coords)} points")
            ax.scatter(coords[:, 0], coords[:, 1], c=colors, s=SCATTERSIZE_VOX, alpha=1)

    for patch in ax.patches[:]:
        new_patch = Polygon(patch.get_xy(), closed=False, edgecolor="grey", facecolor="none", 
        linewidth=LINEWIDTH, zorder=patch.get_zorder() + 2)
        ax.add_patch(new_patch)
        patch.remove()

    # Final row: slice labels
    xlabel = f"CCF {closest_x / 1000:.2f} mm"
    label_ax = fig.add_subplot(gs[1, col])
    label_ax.axis("off")
    label_ax.text(
        0.5, 0.5, 
        xlabel,
        fontsize=12, fontweight="bold", ha="center", va="center", transform=label_ax.transAxes
    )


# --- Colorbars ---
formatter = ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((0, 0))  # always use scientific notation

# MODALITY colorbar vox
cax_prob_vox = fig.add_subplot(gs[0, -1])
cb_prob_vox = fig.colorbar(sm_prob_vox, cax=cax_prob_vox)
cb_prob_vox.ax.tick_params(labelsize=12)  # Tick label fontsize
cb_prob_vox.ax.yaxis.set_tick_params(width=1.5)
cax_prob_vox.set_ylabel(f"{'rel. mass\n' if MODALITY_VOX == 'rel_mass' else 'p'} vox", rotation=90, labelpad=8, fontsize=14, fontweight='bold')
if MODALITY_VOX == "rel_mass":
    cb_prob_vox.formatter = formatter
    cb_prob_vox.update_ticks()
    cb_prob_vox.ax.yaxis.get_offset_text().set_fontweight("bold")
    cb_prob_vox.ax.yaxis.get_offset_text().set_fontsize(12)  
    offset = cb_prob_vox.ax.yaxis.get_offset_text()
    offset.set_x(3.0)  # increase to move further right

    for label in cb_prob_vox.ax.get_yticklabels():
        label.set_fontweight("bold")
else:
    cb_prob_vox.ax.set_yticklabels([t.get_text() for t in cb_prob_vox.ax.get_yticklabels()],
                            fontweight='bold')

# cell prob colorbar
cax_prob_cell = fig.add_subplot(gs[0, -3])
cb_prob_cell = fig.colorbar(sm_prob_cell, cax=cax_prob_cell)
cb_prob_cell.ax.tick_params(labelsize=12)
cb_prob_cell.ax.yaxis.set_tick_params(width=1.5)
cb_prob_cell.ax.set_yticklabels([t.get_text() for t in cb_prob_cell.ax.get_yticklabels()],
                        fontweight='bold')
cax_prob_cell.set_ylabel('p cell', rotation=90, labelpad=8, fontsize=14, fontweight='bold')
cb_prob_cell.ax.yaxis.set_major_formatter(formatter)
cb_prob_cell.update_ticks()


FILENAME_OUT = f"{'_'.join(cluster_descr.split(" "))}_cells_and_vox_{PLAY_FOR_CELLS}um_{MIN_CELL_P_TO_PLOT}p"
plt.tight_layout()
plt.savefig(OUTDIR / f"{FILENAME_OUT}.pdf", dpi=300)
plt.show()
