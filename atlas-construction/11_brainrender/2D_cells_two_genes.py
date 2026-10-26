
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
import json
import brainglobe_heatmap as bgh
from dotenv import dotenv_values

CONFIG = dotenv_values("../../config.env") ## assuming you are in the script on local machine without uv container
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
SCATTERSIZE = params["SCATTERSIZE"] # 0.6 is default
LINEWIDTH = params["LINEWIDTH"] #default 0.6
XLIM = params["XLIM"] # default [-19S00, 1900]
YLIM = params["YLIM"] # default [3200, 800]
GENE_LEFT = params["GENE_LEFT"]
GENE_LEFT_CMAP = params["GENE_LEFT_CMAP"]
GENE_RIGHT = params["GENE_RIGHT"]
GENE_RIGHT_CMAP = params["GENE_RIGHT_CMAP"]
MIN_GE_TO_PLOT_LEFT = params["MIN_GE_TO_PLOT_LEFT"] # exclusive: cells must have GE greater than this
MIN_GE_TO_PLOT_RIGHT = params["MIN_GE_TO_PLOT_RIGHT"] # exclusive: cells must have GE greater than this
RASTERIZE = params.get("RASTERIZE", True)

# derive settings
distances = list(range(START_RANGE, STOP_RANGE, SLICE_THICKNESS))
n_slices = len(distances)
genes = [GENE_LEFT, GENE_RIGHT]

# --- Load Data ---
adata = ad.read_h5ad(SPATIALDIR / "hypomap3d_spatial_cell_level.h5ad", backed="r")
adata = adata[:,genes].to_memory()

df = pd.read_csv(BRAINRENDER_DIR / "coord_space.csv")
df = df.loc[df.z_ccf < 5700,:] # remove voxels from right hemisphere
df.allen_acronym = df.allen_acronym.fillna("HY") # this is a dummy such that regions not in the voxel_store will be included in the heatmap serving as a background

# normalize
sc.pp.normalize_total(adata, target_sum=10_000)
sc.pp.log1p(adata)

# extract data into .obs for convenience
point_array = adata.obsm["coords_ccf_mm_wrap"]
point_array = np.array(point_array, dtype = "float32")
point_array = point_array * 1000 
adata.obs["x_ccf"] = point_array[:, 0]
adata.obs["y_ccf"] = point_array[:, 1]
adata.obs["z_ccf"] = point_array[:, 2]

# filter data to only the slices to be plotted:
adata_mask = np.asarray([False] * len(adata))
df_mask = np.asarray([False] * len(df))
for x in distances:
    closest_x = df.x_ccf.unique()[np.abs(df.x_ccf.unique() - x).argmin()]
    min_x = closest_x - PLAY_FOR_CELLS
    max_x = closest_x + PLAY_FOR_CELLS
    adata_mask |= ((adata.obs.x_ccf > min_x) & (adata.obs.x_ccf < max_x))
    df_mask |= ((df.x_ccf > min_x) & (df.x_ccf < max_x))

adata = adata[adata_mask, :]
df = df.loc[df_mask, :]

# --- Get colormaps while getting vmax from all cells ---
def create_colormap(data, cmap_name):
    vmax = data.max()
    norm = Normalize(vmin=0, vmax=vmax)
    cmap = get_cmap(cmap_name)
    sm = ScalarMappable(norm=norm, cmap=cmap)
    return cmap, norm, sm

gene_colormaps = {}
gene_colormaps[GENE_LEFT] = create_colormap(adata[:, GENE_LEFT].X.toarray().flatten(), GENE_LEFT_CMAP)
gene_colormaps[GENE_RIGHT] = create_colormap(adata[:, GENE_RIGHT].X.toarray().flatten(), GENE_RIGHT_CMAP)

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
cmap_gene_left, norm_gene_left, sm_gene_left = gene_colormaps[GENE_LEFT]
cmap_gene_right, norm_gene_right, sm_gene_right = gene_colormaps[GENE_RIGHT]

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

    gene_slc = adata[(adata.obs.x_ccf > closest_x - PLAY_FOR_CELLS) & (adata.obs.x_ccf < closest_x + PLAY_FOR_CELLS)]
    if gene_slc.shape[0] > 0:
        gene_values_left = gene_slc[:, GENE_LEFT].X.toarray().flatten()
        coords = gene_slc.obs[["z_ccf", "y_ccf"]].values.copy()
        # ALIGN WITH HEATMAP
        _, y_bias, z_bias = scene.slicer.plane0.center_of_mass()
        coords[:, 0] += z_bias
        coords[:, 1] -= y_bias


        # --- RIGHT Hemisphere: (left on the plot) Gene expression ---
        # restrict to MIN_GE_TO_PLOT_LEFT
        min_ge_mask = gene_values_left > MIN_GE_TO_PLOT_LEFT
        gene_values_left = gene_values_left[min_ge_mask]
        coords_ge_left = coords[min_ge_mask, :]

        # sort so highest expression is plotted last
        order = np.argsort(gene_values_left)
        gene_values_left = gene_values_left[order]
        coords_ge_left = coords_ge_left[order]

        colors_left = cmap_gene_left(norm_gene_left(gene_values_left))

        ax.scatter(coords_ge_left[:, 0], coords_ge_left[:, 1], c=colors_left, s=SCATTERSIZE, alpha=1, rasterized=RASTERIZE)

        # --- LEFT Hemisphere: (right on the plot) Gene expression ---
        coords[:, 0] = -coords[:, 0]  # LEFT hemisphere flip
        gene_values_right = gene_slc[:, GENE_RIGHT].X.toarray().flatten()

        # restrict to MIN_GE_TO_PLOT_LEFT
        min_ge_mask = gene_values_right > MIN_GE_TO_PLOT_RIGHT
        gene_values_right = gene_values_right[min_ge_mask]
        coords_ge_right = coords[min_ge_mask, :]

        # sort so highest expression is plotted last
        order = np.argsort(gene_values_right)
        gene_values_right = gene_values_right[order]
        coords_ge_right = coords_ge_right[order]

        colors_right = cmap_gene_right(norm_gene_right(gene_values_right))

        ax.scatter(coords_ge_right[:, 0], coords_ge_right[:, 1], c=colors_right, s=SCATTERSIZE, alpha=1, rasterized=RASTERIZE)

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
# Gene colorbar left
cax_gene_left = fig.add_subplot(gs[0, -3])
cb_gene_left = fig.colorbar(sm_gene_left, cax=cax_gene_left)
cb_gene_left.ax.tick_params(labelsize=12)
cb_gene_left.ax.yaxis.set_tick_params(width=1.5)
cb_gene_left.ax.set_yticklabels([t.get_text() for t in cb_gene_left.ax.get_yticklabels()],
                        fontweight='bold')
cax_gene_left.set_ylabel(GENE_LEFT, rotation=90, labelpad=8, fontsize=14, fontweight='bold')

# Gene colorbar right
cax_gene_right = fig.add_subplot(gs[0, -1])
cb_gene_right = fig.colorbar(sm_gene_right, cax=cax_gene_right)
cb_gene_right.ax.tick_params(labelsize=12)
cb_gene_right.ax.yaxis.set_tick_params(width=1.5)
cb_gene_right.ax.set_yticklabels([t.get_text() for t in cb_gene_right.ax.get_yticklabels()],
                        fontweight='bold')
cax_gene_right.set_ylabel(GENE_RIGHT, rotation=90, labelpad=8, fontsize=14, fontweight='bold')


FILENAME_OUT = f"{START_RANGE}-{STOP_RANGE}_{GENE_LEFT}_{MIN_GE_TO_PLOT_LEFT}_and_{GENE_RIGHT}_{MIN_GE_TO_PLOT_RIGHT}_{PLAY_FOR_CELLS}um"
if RASTERIZE == True:
    FILENAME_OUT += "_rasterized"
plt.tight_layout()
plt.savefig(OUTDIR / f"{FILENAME_OUT}.pdf", dpi=300)
plt.show()
