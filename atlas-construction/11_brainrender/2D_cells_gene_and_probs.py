
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
SCATTERSIZE = params["SCATTERSIZE"] # 0.6 is default
LINEWIDTH = params["LINEWIDTH"] #default 0.6
XLIM = params["XLIM"] # default [-1900, 1900]
YLIM = params["YLIM"] # default [3200, 800]
GENE = params["GENE"]
GENE_CMAP = params["GENE_CMAP"]
PROB_CMAP = params["PROB_CMAP"]
MIN_GE_TO_PLOT = params["MIN_GE_TO_PLOT"] # exclusive: cells must have GE greater than this
MIN_P_TO_PLOT = params["MIN_P_TO_PLOT"] # exclusive: cells must have p greater than this
CLUSTER_TO_PLOT = params["CLUSTER_TO_PLOT"]
MIN_CELL_LOAD_CONFIDENCE = params["MIN_CELL_LOAD_CONFIDENCE"] #float; inclusive ie cells with p greater or equal than this will be loaded
#MIN_P_CELLS_TO_PLOT = params["MIN_P_CELLS_TO_PLOT"] #float; inclusive ie cells with p greater or equal than this will be plotted for cluster probs
RASTERIZE = params.get("RASTERIZE", True)

# derive settings
ID = CLUSTER_TO_PLOT.split(" ")[0] # failsafe for if descr is passed
level_key = f"L{ID[2]}"
distances = list(range(START_RANGE, STOP_RANGE, SLICE_THICKNESS))
n_slices = len(distances)

# --- Load Data ---
adata = ad.read_h5ad(SPATIALDIR / "hypomap3d_spatial_cell_level.h5ad", backed="r")
celltype_idx = list(adata.uns[f"{level_key}_ID"]).index(ID)
load_mask = adata.obsm[f"P_{level_key}"][:, celltype_idx].toarray() >= MIN_CELL_LOAD_CONFIDENCE # load cells of MIN_CELL_CONFIDENCE confidence
adata = adata[load_mask,:].to_memory()

cluster_descr = adata.uns[f"{level_key}_descr"][celltype_idx]

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

adata = adata[adata_mask, :]
df = df.loc[df_mask, :]

# --- Get colormaps while getting vmax from all cells ---
def create_colormap(data, cmap_name):
    vmax = data.max()
    norm = Normalize(vmin=0, vmax=vmax)
    cmap = get_cmap(cmap_name)
    sm = ScalarMappable(norm=norm, cmap=cmap)
    return cmap, norm, sm

prob_colormaps = {}
gene_colormaps = {}

prob_colormaps[ID] = create_colormap(adata.obs[f"{ID}_prob"].values, PROB_CMAP)
if GENE in adata.var_names:
    gene_values = adata[:, GENE].X.toarray().flatten()
    gene_colormaps[ID] = create_colormap(gene_values, GENE_CMAP)

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
# --- Plotting ---
fig = plt.figure(figsize=(n_slices * 4, 2.25))
height_ratios = [1.0, 0.2]
gs = GridSpec(2, n_slices + 3, figure=fig, 
    width_ratios=[1]*n_slices + [0.05, 0.02, 0.05], wspace=0.3,
    hspace=0.01, height_ratios=height_ratios)

# extract colormaps
cmap_prob, norm_prob, sm_prob = prob_colormaps[ID]
cmap_gene, norm_gene, sm_gene = gene_colormaps[ID]

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

    # --- RIGHT Hemisphere: Gene expression ---
    gene_slc = adata[(adata.obs.x_ccf > closest_x - PLAY_FOR_CELLS) & (adata.obs.x_ccf < closest_x + PLAY_FOR_CELLS)]
    if GENE in gene_slc.var_names and gene_slc.shape[0] > 0:
        gene_values = gene_slc[:, GENE].X.toarray().flatten()
        coords = gene_slc.obs[["z_ccf", "y_ccf"]].values.copy()
        # ALIGN WITH HEATMAP
        _, y_bias, z_bias = scene.slicer.plane0.center_of_mass()
        coords[:, 0] += z_bias
        coords[:, 1] -= y_bias

        # restrict to MIN_GE_TO_PLOT
        min_ge_mask = gene_values > MIN_GE_TO_PLOT
        gene_values = gene_values[min_ge_mask]
        coords_ge = coords[min_ge_mask, :]

        # sort so highest expression is plotted last
        order = np.argsort(gene_values)
        gene_values = gene_values[order]
        coords_ge = coords_ge[order]

        colors = cmap_gene(norm_gene(gene_values))

        ax.scatter(coords_ge[:, 0], coords_ge[:, 1], c=colors, s=SCATTERSIZE, alpha=1, rasterized=RASTERIZE)

        # --- LEFT Hemisphere: Cluster probabilities ---
        coords[:, 0] = -coords[:, 0]  # LEFT hemisphere flip

        # restrict to MIN_P_TO_PLOT
        probs_slice = gene_slc.obs[f"{ID}_prob"].values
        min_p_mask = probs_slice > MIN_P_TO_PLOT
        probs_slice = probs_slice[min_p_mask]
        coords_p = coords[min_p_mask]

        # sort so highest probability is plotted last
        order = np.argsort(probs_slice)
        probs_slice = probs_slice[order]
        coords_p = coords_p[order]

        colors = cmap_prob(norm_prob(probs_slice))
        ax.scatter(coords_p[:, 0], coords_p[:, 1], c=colors, s=SCATTERSIZE, alpha=1, rasterized=RASTERIZE)

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
# Probability colorbar
cax_prob = fig.add_subplot(gs[0, -1])
cb_prob = fig.colorbar(sm_prob, cax=cax_prob)
cb_prob.ax.tick_params(labelsize=12)  # Tick label fontsize
cb_prob.ax.yaxis.set_tick_params(width=1.5)
cb_prob.ax.set_yticklabels([t.get_text() for t in cb_prob.ax.get_yticklabels()],
                        fontweight='bold')
cax_prob.set_ylabel(f"prob.", rotation=90, labelpad=8, fontsize=14, fontweight='bold')

# Gene colorbar
cax_gene = fig.add_subplot(gs[0, -3])
cb_gene = fig.colorbar(sm_gene, cax=cax_gene)
cb_gene.ax.tick_params(labelsize=12)
cb_gene.ax.yaxis.set_tick_params(width=1.5)
cb_gene.ax.set_yticklabels([t.get_text() for t in cb_gene.ax.get_yticklabels()],
                        fontweight='bold')
cax_gene.set_ylabel(GENE, rotation=90, labelpad=8, fontsize=14, fontweight='bold')

FILENAME_OUT = f"{'_'.join(cluster_descr.split(" "))}_and_{GENE}_{PLAY_FOR_CELLS}um_{MIN_P_TO_PLOT}p_{MIN_GE_TO_PLOT}ge"
if RASTERIZE == True:
    FILENAME_OUT += "_rasterized"
plt.tight_layout()
plt.savefig(OUTDIR / f"{FILENAME_OUT}.pdf", dpi=300)
plt.show()
