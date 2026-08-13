
import os 
import argparse
from pathlib import Path
import json
import zarr
import numpy as np
import pandas as pd
from brainrender import Scene, settings
from brainrender.actors import Volume
import vedo
from vedo import settings as vsettings
from matplotlib.colors import LinearSegmentedColormap, Normalize
from vedo import build_lut
from vedo import Volume as vVolume
from vedo.addons import ScalarBar
from skimage.filters import threshold_otsu
from dotenv import dotenv_values

import brainrender
import vtk 

print(brainrender.__version__)
print(vedo.__version__)
print(vtk.vtkVersion().GetVTKVersion())

settings.ROOT_ALPHA = 0.1

CONFIG = dotenv_values("../../config.env") # assuming you are in the script on local machine without uv container
DATADIR = Path(CONFIG['DATADIR'])

SPATIALDIR = DATADIR / "9_objects" / "spatial"
BRAINRENDER_DIR = DATADIR / "11_brainrender" 
OUTDIR = BRAINRENDER_DIR / "3D_plots"
os.makedirs(OUTDIR, exist_ok=True)

ZARR_FP = SPATIALDIR / "hypomap3d_spatial_voxel_level.zarr"

# brainrender params
X_RANGE = (0, 13200)
Y_RANGE = (0, 8000)
Z_RANGE = (0, 11400)
MAX_Y = 4000 # to cut off SFO from HY mesh
VOXEL_SIZE = 25  # in microns

parser = argparse.ArgumentParser()
parser.add_argument(
    "params",
    help="Path to parameters JSON file."
)

args = parser.parse_args()

with open(args.params) as f:
    params = json.load(f)

CAMERA_POS = params["camera_pos"]
SHOW_ROOT = params["show_root"]

RESTRICT_P_TO_HY = params["restrict_p_to_hy"]
USE_OTSU = params["use_otsu"]
USE_MIN_P = params["use_min_p"]

if USE_MIN_P:
    min_p = config["min_p"]

MODALITY = params["modality"]
CLUSTERS_TO_PLOT = params["clusters_to_plot"]
voxel_store = zarr.open(SPATIALDIR / 'hypomap3d_spatial_voxel_level.zarr')


def get_camera(view="sagittal", show_root=True):
    if show_root:
        center = (6600, 4000, -5700)
        if view == "coronal":
            return dict(pos=(-28000, 4000, -5700), focalPoint=center, viewup=(0, -1, 0))
        elif view == "sagittal":
            return dict(pos=(6600, 4000, 33000), focalPoint=center, viewup=(0, -1, 0))
        elif view == "axial":
            return dict(pos=(6600, -45000, -5700), focalPoint=center, viewup=(0, 0, -1))
        if view == "frontleft":
            return dict(pos=(-17000, 4000, 25000), focalPoint=center, viewup=(0, -1, 0))
        
    else:
        center = (6600, 6000, -5700)
        if view == "coronal":
            return dict(pos=(-3000, 6000, -5700), focalPoint=center, viewup=(0, -1, 0))
        elif view == "sagittal":
            return dict(pos=(6600, 6000, 10000), focalPoint=center, viewup=(0, -1, 0))
        elif view == "axial":
            #center = (6600, 6000, -4700)
            return dict(pos=(6600, -6000, -5700), focalPoint=center, viewup=(0, 0, -1))
        if view == "frontleft":
            return dict(pos=(-3000, 6000, 1000), focalPoint=center, viewup=(0, -1, 0))


df = pd.read_csv(BRAINRENDER_DIR / "coord_space.csv")

for CLUSTER_TO_PLOT in CLUSTERS_TO_PLOT:

    print(f"Processing {CLUSTER_TO_PLOT}")

    for camera_pos in CAMERA_POS: # you need to create new scenes for each camera position
        print(f"Getting camera position: {camera_pos}")

        # Load data to be plotted
        ID = CLUSTER_TO_PLOT.split(" ")[0]
        level_key = f"L{ID[2]}"
        idx = voxel_store[level_key]["ids"][:].tolist().index(ID)
        probs = voxel_store[level_key][MODALITY][idx, :, :, :].reshape(-1)
        probs = np.nan_to_num(probs, nan=0.0)
        cluster_descr = voxel_store[level_key]["descr"][:].tolist()[idx]
        df["probability"] = probs
        flag = voxel_store[level_key]["flags"][:].tolist()[idx]
        print(f"Flag: {flag}")

        if RESTRICT_P_TO_HY: # set all probabilities outside the hy to 0 
            mask = df["is_inside_HY"].values  
            df.loc[~mask, "probability"] = 0

        if df.probability.sum() == 0:
            print(f"Skipping {CLUSTER_TO_PLOT} for being zero everywhere")
            continue

        print(f"Creating voxel grid")
        # Step 2: Compute shape of full volume
        shape_x = (X_RANGE[1] - X_RANGE[0]) // VOXEL_SIZE
        shape_y = (Y_RANGE[1] - Y_RANGE[0]) // VOXEL_SIZE
        shape_z = (Z_RANGE[1] - Z_RANGE[0]) // VOXEL_SIZE
        volume_shape = (shape_x, shape_y, shape_z)

        # Step 3: Initialize grid with zeros (padding)
        grid = np.full(volume_shape, fill_value=0, dtype=np.float32)

        # Step 4: Convert CCF coordinates to voxel indices
        x_idx = ((df['z_ccf'] - X_RANGE[0]) / VOXEL_SIZE).astype(int) # notice z_ccf here. this is correct if as_surface=False, but incorrect if as_surface=True
        y_idx = ((df['y_ccf'] - Y_RANGE[0]) / VOXEL_SIZE).astype(int)
        z_idx = ((df['x_ccf'] - Z_RANGE[0]) / VOXEL_SIZE).astype(int)

        # Step 5: Fill in probabilities
        grid[x_idx, y_idx, z_idx] = df['probability'].values

        print("Creating scene")
        scene = Scene(
            root=SHOW_ROOT
        )

        hy = scene.add_brain_region('HY', alpha=0.1, color="red", hemisphere="right"
            )

        cut_plane = scene.atlas.get_plane(pos=(0, MAX_Y, 0), norm=(0, 1, 0))

        scene.slice(cut_plane, actors=[hy])

        if not USE_OTSU:
            if not USE_MIN_P:
                if RESTRICT_P_TO_HY:
                    min_perc = 96
                    #min_perc = 95
                else:
                    min_perc = 98 # use this percentile as the min p value in the plot
                try:
                    min_p = np.percentile(grid[grid > 0], min_perc)
                except:
                    min_p = 0
                print(f"{min_p=}")

        if not USE_OTSU:

            if not USE_MIN_P:

                lut_transparent = build_lut(
                    [
                        #(0, "white", 0.0),
                        #(np.percentile(grid[grid > 0], 97.99999), "white", 0.00),
                        (np.percentile(grid[grid > 0], min_perc), "blue", 0.001), 
                        (np.percentile(grid[grid > 0], 99.9), "green", 0.01),
                        (np.percentile(grid[grid > 0], 99.95), "yellow", 0.6),
                        (np.max(grid), "white", 0.99),
                    ],
                    interpolate=True,
            )

            else:
                #print("get lut")
                lut_transparent = build_lut(
                [
                    #(0, "white", 0.0),
                    #(min_p - 0.000001, "white", 0.00),
                    (min_p, "blue", 0.001),
                    (np.percentile(grid[grid > 0], 99.9), "green", 0.01),
                    (np.percentile(grid[grid > 0], 99.95), "yellow", 0.6),
                    (np.max(grid), "white", 0.99),
                ],
                interpolate=True,
            )

        else:
            min_p = threshold_otsu(grid)

            lut_transparent = build_lut(
                [
                    #(0, "white", 0.0),
                    #(min_p - 0.000001, "white", 0.00),
                    (min_p, "blue", 0.001),
                    (np.percentile(grid[grid > 0], 99.9), "green", 0.01),
                    (np.percentile(grid[grid > 0], 99.95), "yellow", 0.6),
                    (np.max(grid), "white", 0.99),
                ],
                interpolate=True,
            )


        # set all values lower than min_p to min_p. according to color scale these values will all be transparent, so not displayed.
        #print("test")
        grid[grid < min_p] = min_p

        print("Adding Volume actor")
        actor = vVolume(
            grid,
            spacing=[25, 25, 25],
        )

        actor.permute_axes(2, 1, 0)
        mtx = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, 1]]
        actor.apply_transform(mtx)


        actor.cmap(lut_transparent, vmin = min_p # dont use min_p in cases when your LUT starts with white up to min_p. this seems to let the white part start at min_p instead 
        ) # overwrite lut

        print("Adding scalar bar")
        sc2d = ScalarBar(
            actor,
            # title="Membership probability",
            c="k",
            size=[0, 1],
            title_yoffset=20,
            font_size=20,
            #label_size=1.8,
            use_alpha = False,
            nlabels=5,
            label_format = "%.3f" if MODALITY == "prob" else f"%.1e"
        )
        # colorscale formatting
        # pos
        sc2d.SetPosition(0.8, 0.05)   # center-bottom
        sc2d.SetWidth(0.13)
        sc2d.SetHeight(0.88)
        sc2d.SetTextPositionToSucceedScalarBar()   
        # labels
        label_prop = sc2d.GetLabelTextProperty()
        label_prop.SetFontSize(80 if MODALITY == "prob" else 55)
        label_prop.SetBold(True)


        # Extract brainrender actors
        brainrender_actors = scene.actors  # This is a list of vedo actors
        all_actors = brainrender_actors + [sc2d, actor] 


        # get camera
        hy_com = hy.center_of_mass()
        hy_com[2] *= -1

        pos = hy_com.copy()
        pos[0] -= 7000

        if camera_pos == "frontleft":
            pos[2] += 4500 

        if camera_pos == "coronal":
            hy_com[2] += 1200

        if camera_pos == "sagittal":
            pos[0] += hy_com[0]
            pos[2] -= hy_com[2] * 2

        camera = {
            "pos": pos,     # camera position in world coordinates
            "focalPoint": hy_com,     # point camera looks at
            "viewup": (0, -1, 0),         # "up" direction for the camera: -y ccf
        }


        screenshot_path = str(OUTDIR) + f"/{'_'.join(cluster_descr.split(' '))}_3D_{camera_pos}_{MODALITY}"
        if SHOW_ROOT:
            screenshot_path += "_with_root"

        if RESTRICT_P_TO_HY:
            screenshot_path += "_hy"

        if USE_OTSU:
            screenshot_path += "_otsu"

        if USE_MIN_P:
            screenshot_path += f"_minp{min_p}"

        screenshot_path += ".svg"

        vedo.show(*all_actors, size=(2000, 2000), camera=camera,
        screenshot=screenshot_path,
        )

        del scene

        print(f"Saved at {screenshot_path}")