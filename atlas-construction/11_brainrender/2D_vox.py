from pathlib import Path
import argparse
import zarr
import os 
import pandas as pd 
import brainglobe_heatmap as bgh
import brainrender
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.colors import Normalize, to_hex
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Polygon
from matplotlib.gridspec import GridSpec
import json
import gc
from dotenv import dotenv_values

# Set cartoon shader
brainrender.settings.SHADER_STYLE = "cartoon"

CONFIG = dotenv_values("../../config.env") # assuming you are in the script on local machine without uv container
DATADIR = Path(CONFIG['DATADIR'])

SPATIALDIR = DATADIR / "9_objects" / "spatial"
BRAINRENDER_DIR = DATADIR / "11_brainrender" 
OUTDIR = BRAINRENDER_DIR / "2D_plots"
os.makedirs(OUTDIR, exist_ok=True)

ZARR_FP = SPATIALDIR / "hypomap3d_spatial_voxel_level.zarr"


parser = argparse.ArgumentParser()
parser.add_argument(
    "params",
    help="Path to parameters JSON file."
)

args = parser.parse_args()

with open(args.params) as f:
    params = json.load(f)
    

# Parameters
SLICE_THICKNESS = params["SLICE_THICKNESS"]
START_RANGE = params["START_RANGE"] # 4750 is default, 6200 for periARH, 6500 for tany
STOP_RANGE = params["STOP_RANGE"] # 8500 is default, 7800 for periARH, 7800 for tany
ANNOTATE_REGIONS = params["ANNOTATE_REGIONS"] # False is default
SCATTERSIZE = params["SCATTERSIZE"] # 0.6 is default
LINEWIDTH = params["LINEWIDTH"] #default 0.2
XLIM = params["XLIM"] # default [-100, 1900]
YLIM = params["YLIM"] # default [3200, 800]
MIN_P = params["MIN_P"] # default 0.001
#MIN_P = 0
CMAP_NAME = params["CMAP_NAME"] # default Blues
CLUSTERS_TO_PLOT = params["CLUSTERS_TO_PLOT"]
RESTRICT_P_TO_HY = params["RESTRICT_P_TO_HY"]
MODALITY = params["MODALITY"]
IS_BREGMA_COORDS = params["IS_BREGMA_COORDS"]

MIN_P_THR_QUANTILE = False # whether to use percentiles as min_p
PLOT_NUCLEI_BELOW = False

ALPHA = 0.2 if PLOT_NUCLEI_BELOW else 1
FACECOLOR = "white" if PLOT_NUCLEI_BELOW else "none"

SCALE_FACTOR = params.get("SCALE_FACTOR", 0.4)

# Helper for heatmap values
def get_values_for_slice(slice_df):
    return {k: 0 for k in slice_df.allen_acronym.unique()}
# helper for Bregma --> CCF conversion 
def precompute_ccf_bregma_matches(voxel_grid):
    
    n = voxel_grid.shape[0] # slice along coronal (x)
    
    bound_voxel_coord = np.arange(n)
    
    voxel_coord = bound_voxel_coord + 175 # FUTURE: not hardcoded, from bounding box instead
    
    ccf_coord = voxel_coord * 25  # um
    
    midpoint_y = (voxel_grid.shape[1] + 158) / 2 # FUTURE: derive 158 from bounding box object
    midpoint_y = midpoint_y * 25  # um
    
    theta = np.deg2rad(5)

    bregma = (5400, 440, 5700)
    
    adj_x =  (ccf_coord - bregma[0]) * -1
    adj_y =  (midpoint_y - bregma[1]) * 0.9434
    
    bregma_coord = (adj_x * np.cos(theta)) - (adj_y * np.sin(theta))
    
    return np.column_stack((ccf_coord, bregma_coord))


# Load data
df = pd.read_csv(BRAINRENDER_DIR / "coord_space.csv")
#df = df.fillna("0")
na_mask = df.allen_id == 0
df = df.loc[~na_mask, :]
df.allen_acronym = df.allen_acronym.fillna("HY") # this is a dummy such that regions not in the voxel_store will be included in the heatmap serving as a background
voxel_store = zarr.open(SPATIALDIR / 'hypomap3d_spatial_voxel_level.zarr')

if IS_BREGMA_COORDS:
    matches = precompute_ccf_bregma_matches(voxel_store["region"]["region_labels"])
    ccf_coords = matches[:, 0]
    bregma_coords = matches[:, 1]

for CLUSTER_TO_PLOT in CLUSTERS_TO_PLOT:

    print(f"Processing {CLUSTER_TO_PLOT}")

    # Load data to be plotted
    ID = CLUSTER_TO_PLOT.split(" ")[0]
    level_key = f"L{ID[2]}"
    idx = voxel_store[level_key]["ids"][:].tolist().index(ID)
    probs = voxel_store[level_key][MODALITY][idx, :, :, :].reshape(-1)
    probs = np.nan_to_num(probs, nan=0.0)
    probs = probs[~na_mask]
    cluster_descr = voxel_store[level_key]["descr"][:].tolist()[idx]
    df["probability"] = probs
    flag = voxel_store[level_key]["flags"][:].tolist()[idx]
    print(f"Flag: {flag}")

    if RESTRICT_P_TO_HY:
        mask = df["is_inside_HY"].values 
        df.loc[~mask, "probability"] = 0

    distances = np.arange(
        START_RANGE,
        STOP_RANGE,
        SLICE_THICKNESS
    ).tolist()
    if IS_BREGMA_COORDS:
        distances = distances[::-1]
    n_slices = len(distances)

    if np.sum(df["probability"].values) == 0:
        print(f"Skipping {cluster_descr}. All probs are zero")
        continue

    fig = plt.figure(figsize=(n_slices * 2 * SCALE_FACTOR, 3.5 * SCALE_FACTOR))  # One row per cluster plot
    height_ratios = [0.2, 1.0, 0.2]
    gs = GridSpec(3, n_slices, figure=fig, wspace=0.1, hspace=0.01, height_ratios=height_ratios)

    # calculate a min probability to be plotted based on a quantile of the nonzero probabilities
    if MIN_P_THR_QUANTILE:
        if RESTRICT_P_TO_HY:
            min_perc = 96
        else:
            min_perc = 98 # use this percentile as the min p value in the plot
        MIN_P = np.percentile(df.loc[df["probability"]>0, "probability"].values, min_perc)

    # Cluster title
    title_ax = fig.add_subplot(gs[0, :])
    title_ax.axis("off")
    titletext = f"{cluster_descr}"
    title_ax.text(
        0.5, 0.5, titletext,
        fontsize=15 * SCALE_FACTOR, ha="center", va="center", weight="bold", transform=title_ax.transAxes
    )

    closest_xs = []
    unique_x = df.x_ccf.unique()

    for j, distance in enumerate(distances):
        if not IS_BREGMA_COORDS:
            closest_x = df.x_ccf.unique()[np.abs(df.x_ccf.unique() - distance).argmin()]
        else:
            idx = np.abs(bregma_coords - distance).argmin()
            ccf = ccf_coords[idx]
            closest_x = unique_x[np.abs(unique_x - ccf).argmin()]

        closest_xs.append(closest_x) # keep for color map
        spatial_slice = df.loc[(df.x_ccf == closest_x), :]

        scene = bgh.Heatmap(
            get_values_for_slice(spatial_slice),
            position=closest_x,
            orientation="frontal",
            thickness=1,
            format="2D",
            annotate_regions=ANNOTATE_REGIONS,
            title=f"X_CCF = {closest_x/1000} mm",
            hemisphere="both",
            root=False,
        )

        # allen_color_mapper = spatial_slice.set_index("allen_acronym")["allen_color"].to_dict()
        # for k in allen_color_mapper:
        #     scene.colors[k] = "#FFFFFF"

        _, y_bias, z_bias = scene.slicer.plane0.center_of_mass()

        ax = fig.add_subplot(gs[1, j])
        scene.plot_subplot(fig=fig, ax=ax, show_cbar=False, hide_axes=True)

        if PLOT_NUCLEI_BELOW:
            for patch in ax.patches[:]:
                new_patch = Polygon(
                    patch.get_xy(), closed=False,
                    edgecolor="grey", facecolor=FACECOLOR,
                    linewidth=LINEWIDTH, zorder=patch.get_zorder()
                )
                ax.add_patch(new_patch)
                patch.remove()

        cluster_slice = df.loc[
            (df.x_ccf == closest_x) & 
            (df.z_ccf < 5700) & 
            (df["probability"] >= MIN_P),
            ["x_ccf", "y_ccf", "z_ccf", "probability"]
        ]

        probs = cluster_slice["probability"].values
        coords_2D = cluster_slice[["z_ccf", "y_ccf"]].values

        if not len(probs) == 0:
            coords_2D[:, 0] += z_bias
            coords_2D[:, 1] -= y_bias
            coords_2D[:, 0] = -coords_2D[:, 0]

            norm = Normalize(vmin=MIN_P, vmax=df.loc[df.x_ccf.isin(closest_xs), "probability"].max())
            cmap = plt.get_cmap(CMAP_NAME)
            sm = ScalarMappable(norm=norm, cmap=cmap)
            colors = [to_hex(cmap(norm(p))) for p in probs]

            scatter = ax.scatter(coords_2D[:, 0], coords_2D[:, 1], c=colors, s=SCATTERSIZE * SCALE_FACTOR, alpha=ALPHA)
            scatter.set_rasterized(True)

        if not PLOT_NUCLEI_BELOW:
            for patch in ax.patches[:]:
                new_patch = Polygon(
                    patch.get_xy(), closed=False,
                    edgecolor="grey", 
                    #edgecolor="white", 
                    facecolor=FACECOLOR,
                    linewidth=LINEWIDTH * SCALE_FACTOR, zorder=patch.get_zorder() + 2
                )
                ax.add_patch(new_patch)
                patch.remove()

        ax.set_rasterized(True)
        ax.set_xlim(XLIM)
        ax.set_ylim(YLIM)
        ax.set_title("")
        ax.set_ylabel("")

        if j == n_slices - 1:
            cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
            cbar.ax.tick_params(labelsize=15 * SCALE_FACTOR, labelcolor='black')
            for label in cbar.ax.get_yticklabels():
                label.set_fontweight('bold') 
            if MODALITY == "rel_mass":
                from matplotlib.ticker import ScalarFormatter
                formatter = ScalarFormatter(useMathText=True)
                formatter.set_scientific(True)
                formatter.set_powerlimits((0, 0))  # always use scientific notation
                cbar.formatter = formatter
                cbar.update_ticks()
                cbar.ax.yaxis.get_offset_text().set_fontweight("bold")
                cbar.ax.yaxis.get_offset_text().set_fontsize(12 * SCALE_FACTOR)  
                offset = cbar.ax.yaxis.get_offset_text()
                offset.set_x(8.0)  # increase to move further right
                cbar.ax.set_title(
                "rel. m",
                fontsize=15 * SCALE_FACTOR,
                fontweight="bold",
                pad=4,
                x=5,   # shift to the right
            )
            else:
                cbar.ax.set_title("p", fontsize=15 * SCALE_FACTOR, fontweight='bold', pad=4)

        del scene # debug


    # Final row: slice labels
    for j, distance in enumerate(distances):
        if not IS_BREGMA_COORDS:
            closest_x = df.x_ccf.unique()[np.abs(df.x_ccf.unique() - distance).argmin()]
            xlabel = f"CCF {closest_x / 1000:.2f} mm"
            fontsize=12 * SCALE_FACTOR
        else:
            idx = np.abs(bregma_coords - distance).argmin()
            ccf = ccf_coords[idx]
            closest_x = unique_x[np.abs(unique_x - ccf).argmin()]
            xlabel = f"CCF {closest_x / 1000:.2f} mm\nBregma ~ {distance / 1000:.2f}"
            fontsize=10 * SCALE_FACTOR
        label_ax = fig.add_subplot(gs[2, j])
        label_ax.axis("off")
        label_ax.text(
            0.5, 0.5, 
            xlabel,
            fontsize=fontsize, fontweight="bold", ha="center", va="center", transform=label_ax.transAxes
        )


    plt.tight_layout()

    # Save with cluster name suffix
    fig_save_path = str(OUTDIR) + f"/{'_'.join(cluster_descr.split(' '))}_2D_{MODALITY}"
    if RESTRICT_P_TO_HY:
        fig_save_path += "_hyP"
    if MIN_P_THR_QUANTILE:
        fig_save_path += "_percentile"
    fig.savefig(fig_save_path + f"_{SLICE_THICKNESS}.pdf", dpi=300)
    plt.close(fig)

    gc.collect()
