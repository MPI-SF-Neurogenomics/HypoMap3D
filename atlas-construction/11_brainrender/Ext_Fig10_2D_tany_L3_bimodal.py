
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
import zarr
import brainglobe_heatmap as bgh
from dotenv import dotenv_values

CONFIG = dotenv_values("../../config.env") # assuming you are in the script on local machine without uv container
DATADIR = Path(CONFIG['DATADIR'])

SPATIALDIR = DATADIR / "9_objects" / "spatial"
BRAINRENDER_DIR = DATADIR / "11_brainrender" 
OUTDIR = BRAINRENDER_DIR / "2D_plots"
os.makedirs(OUTDIR, exist_ok=True)

ZARR_FP = SPATIALDIR / "hypomap3d_spatial_voxel_level.zarr"

MODALITY = "rel_mass"
FILENAME_OUT = f"tany_{MODALITY}_and_genes"

slice_thickness = 300 #100
play_for_cells = slice_thickness / 2 
xlim = (-500, 500)
ylim = (3200, 2000)
scattersize_probs = 30
scattersize_cells = 6
linewidth = 0.6
annotate_regions = False

distances_start = 6500
distances_stop = 7800
distances = list(range(distances_start, distances_stop, slice_thickness))

# --- Define Clusters & Genes ---
cluster_gene_map = {
    "NN3-27 Tanycyte_alpha1": "Mafb",
    "NN3-26 Tanycyte_alpha2": "Vcan", #"Pdzph1", #"Cadm2",
    "NN3-25 Tanycyte_beta1": "Frzb",
    "NN3-24 Tanycyte_beta2": "Scn7a"
}

cluster_palettes_map = {
    "NN3-27 Tanycyte_alpha1": "Blues",
    "NN3-26 Tanycyte_alpha2": "Reds", 
    "NN3-25 Tanycyte_beta1": "Purples",
    "NN3-24 Tanycyte_beta2": "Greens"
}

# --- Load Data ---
adata = ad.read_h5ad(SPATIALDIR / "hypomap3d_spatial_cell_level.h5ad", backed="r")
tany_idx = list(adata.uns["L2_descr"]).index("NN2-13 Tanycyte")
load_mask = adata.obsm["P_L2"][:, tany_idx].toarray() > 0.5 # high conf tanys
adata = adata[load_mask,:].to_memory()

voxel_store = zarr.open(SPATIALDIR / 'hypomap3d_spatial_voxel_level.zarr')

df = pd.read_csv(BRAINRENDER_DIR / "coord_space.csv")

clusters = list(cluster_gene_map.keys())
for cl in clusters:

    idx = voxel_store["L3"]["descr"][:].tolist().index(cl)
    probs = voxel_store["L3"][MODALITY][idx, :, :, :].reshape(-1)
    probs = np.nan_to_num(probs, nan=0.0)
    df[cl] = probs

df = df.loc[df.z_ccf < 5700,:] # remove voxels from right hemisphere
df.allen_acronym = df.allen_acronym.fillna("HY") # this is a dummy such that regions not in the voxel_store will be included in the heatmap serving as a background

# normalize
sc.pp.normalize_total(adata, target_sum=10_000)
sc.pp.log1p(adata)

point_array = adata.obsm["coords_ccf_mm_wrap"]
point_array = np.array(point_array, dtype = "float32")
point_array = point_array * 1000 
adata.obs["x_ccf"] = point_array[:, 0]
adata.obs["y_ccf"] = point_array[:, 1]
adata.obs["z_ccf"] = point_array[:, 2]

# --- Colormaps ---
def create_colormap(data, cmap_name):
    vmax = data.max()
    norm = Normalize(vmin=0, vmax=vmax)
    cmap = get_cmap(cmap_name)
    sm = ScalarMappable(norm=norm, cmap=cmap)
    return cmap, norm, sm


prob_colormaps = {}
gene_colormaps = {}

for cluster in clusters:
    prob_colormaps[cluster] = create_colormap(df[cluster], cluster_palettes_map[cluster])

    gene = cluster_gene_map[cluster]
    if gene in adata.var_names:
        gene_values = adata[:, gene].X.toarray().flatten()
        gene_colormaps[cluster] = create_colormap(gene_values, "viridis")

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
        annotate_regions=annotate_regions,
        hemisphere="both",
        root=False,
    )

    scenes_by_x[closest_x] = scene

# --- Plotting ---
fig = plt.figure(figsize=(len(distances) * 4, len(clusters) * 4))
gs = GridSpec(len(clusters), len(distances) + 3, figure=fig, width_ratios=[1]*len(distances) + [0.05, 0.01, 0.05], wspace=0.3)

for row, cluster in enumerate(clusters):
    gene = cluster_gene_map[cluster]
    cmap_prob, norm_prob, sm_prob = prob_colormaps[cluster]
    cmap_gene, norm_gene, sm_gene = gene_colormaps[cluster]

    if MODALITY == "rel_mass":
        highpass = np.percentile(df[cluster].values, 99)

    for col, x in enumerate(distances):
        closest_x = df.x_ccf.unique()[np.abs(df.x_ccf.unique() - x).argmin()]
        scene = scenes_by_x[closest_x]
        ax = fig.add_subplot(gs[row, col])
        scene.plot_subplot(fig=fig, ax=ax, show_cbar=False, hide_axes=True)

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect('equal')

        # --- LEFT Hemisphere: Cluster probabilities ---
        if MODALITY == "prob":
            slc = df[(df.x_ccf == closest_x) & (df[cluster] > 0.001)]
        else:
            slc = df[(df.x_ccf == closest_x) & (df[cluster] > highpass)]
        if not slc.empty:
            coords = slc[["z_ccf", "y_ccf"]].values.copy()
            prob_values = slc[cluster].values
            # sort so highest expression is plotted last
            order = np.argsort(prob_values)
            prob_values = prob_values[order]
            coords = coords[order]
            colors = cmap_prob(norm_prob(prob_values))

            # ALIGN WITH HEATMAP
            _, y_bias, z_bias = scene.slicer.plane0.center_of_mass()
            coords[:, 0] += z_bias
            coords[:, 1] -= y_bias
            coords[:, 0] = -coords[:, 0]  # LEFT hemisphere flip

            #print(f"[{cluster}] LEFT @ x={closest_x}: {len(coords)} points")
            ax.scatter(coords[:, 0], coords[:, 1], c=colors, s=scattersize_probs, alpha=1)

        # --- RIGHT Hemisphere: Gene expression ---
        gene_slc = adata[(adata.obs.x_ccf > closest_x - play_for_cells) & (adata.obs.x_ccf < closest_x + play_for_cells)]
        if gene in gene_slc.var_names and gene_slc.shape[0] > 0:
            gene_values = gene_slc[:, gene].X.toarray().flatten()
            coords = gene_slc.obs[["z_ccf", "y_ccf"]].values.copy()
            colors = cmap_gene(norm_gene(gene_values))
            # ALIGN WITH HEATMAP
            _, y_bias, z_bias = scene.slicer.plane0.center_of_mass()
            coords[:, 0] += z_bias
            coords[:, 1] -= y_bias
            # RIGHT hemisphere: no flip

            #print(f"[{cluster}] RIGHT @ x={closest_x}: {len(coords)} points")
            ax.scatter(coords[:, 0], coords[:, 1], c=colors, s=scattersize_cells, alpha=1)

        for patch in ax.patches[:]:
            new_patch = Polygon(patch.get_xy(), closed=False, edgecolor="grey", facecolor="none", linewidth=linewidth, zorder=patch.get_zorder() + 2)
            ax.add_patch(new_patch)
            patch.remove()

    # --- Colorbars ---
    # Probability colorbar
    cax_prob = fig.add_subplot(gs[row, -1])
    cb_prob = fig.colorbar(sm_prob, cax=cax_prob)
    cb_prob.ax.tick_params(labelsize=20)  # Tick label fontsize
    cb_prob.ax.yaxis.set_tick_params(width=1.5)
    cb_prob.ax.set_yticklabels([t.get_text() for t in cb_prob.ax.get_yticklabels()],
                            fontweight='bold')
    if MODALITY == "rel_mass":
        from matplotlib.ticker import ScalarFormatter
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((0, 0))  # always use scientific notation
        cb_prob.formatter = formatter
        cb_prob.update_ticks()
        cb_prob.ax.yaxis.get_offset_text().set_fontweight("bold")
        cb_prob.ax.yaxis.get_offset_text().set_fontsize(20)  
        offset = cb_prob.ax.yaxis.get_offset_text()
        offset.set_x(5.0)  # increase to move further right
        cax_prob.set_ylabel(f"rel. m {cluster.split("_")[-1]}", rotation=90, labelpad=15, fontsize=22, fontweight='bold')
    else:
        cax_prob.set_ylabel(f"p {cluster.split("_")[-1]}", rotation=90, labelpad=15, fontsize=22, fontweight='bold')

    # Gene colorbar
    cax_gene = fig.add_subplot(gs[row, -3])
    cb_gene = fig.colorbar(sm_gene, cax=cax_gene)
    cb_gene.ax.tick_params(labelsize=20)
    cb_gene.ax.yaxis.set_tick_params(width=1.5)
    cb_gene.ax.set_yticklabels([t.get_text() for t in cb_gene.ax.get_yticklabels()],
                            fontweight='bold')
    cax_gene.set_ylabel(gene, rotation=90, labelpad=15, fontsize=22, fontweight='bold')

plt.tight_layout()
plt.savefig(OUTDIR / f"{FILENAME_OUT}.pdf", dpi=300)
plt.show()
