
import os 
import argparse
import json 
from pathlib import Path
import numpy as np
import pandas as pd
import brainrender
from brainrender import Scene
from brainrender.actors import Points, PointsDensity, Volume
from vedo import Plotter  
from vedo.colors import color_map
from vedo import ScalarBar, show
from vedo import Points as vPoints
from dotenv import dotenv_values
from vedo.addons import ScalarBar
import anndata as ad 
from dotenv import dotenv_values

# Set cartoon shader
# #brainrender.settings.SHADER_STYLE = "off" #"plastic" #[cartoon, metallic, plastic, shiny, glossy]
# brainrender.settings.DEFAULT_CAMERA = "frontal" #"three_quarters" #"frontal" if USE_VEDO = False
CONFIG = dotenv_values("../../config.env") # assuming you are in the script on local machine without uv container
DATADIR = Path(CONFIG['DATADIR'])

SPATIALDIR = DATADIR / "9_objects" / "spatial"
BRAINRENDER_DIR = DATADIR / "11_brainrender" 
OUTDIR = BRAINRENDER_DIR / "3D_plots"
os.makedirs(OUTDIR, exist_ok=True)

MAX_Y = 4000 
USE_VEDO = False

parser = argparse.ArgumentParser()
parser.add_argument(
    "params",
    help="Path to parameters JSON file."
)

args = parser.parse_args()

with open(args.params) as f:
    params = json.load(f)

CAMERA_POS = params["camera_pos"] #list
SHOW_ROOT = params["show_root"]
RESTRICT_P_TO_HY = params["restrict_p_to_hy"]
MIN_P = params["min_p"] 
CLUSTERS_TO_PLOT = params["clusters_to_plot"] #list

adata = ad.read_h5ad(SPATIALDIR / "hypomap3d_spatial_cell_level.h5ad", backed="r")

point_array = adata.obsm["coords_ccf_mm_wrap"]
point_array = np.array(point_array, dtype = "float32")
point_array = point_array * 1000 

for cluster_to_plot in CLUSTERS_TO_PLOT:

    for camera_pos in CAMERA_POS:

        level_key = f"L{cluster_to_plot[2]}"
        cluster_id = cluster_to_plot.split(" ")[0]
        idx = list(adata.uns[f"{level_key}_ID"]).index(cluster_id)
        cluster_descr = adata.uns[f"{level_key}_descr"][idx]
        probs = adata.obsm[f"P_{level_key}"][:, idx].toarray().ravel()

        mask = probs > MIN_P 
        cells = point_array[mask] 
        probs = probs[mask]

        scene = Scene(atlas_name='allen_mouse_25um', root=False)
        if SHOW_ROOT:
            scene.add_brain_region('root', alpha=0.1, hemisphere="right")

        hy = scene.add_brain_region('HY', alpha=0.1, color="red", hemisphere="right"
            )

        cut_plane = scene.atlas.get_plane(pos=(0, MAX_Y, 0), norm=(0, 1, 0))

        scene.slice(cut_plane, actors=[hy])

        if RESTRICT_P_TO_HY:
            ids = hy.mesh.inside_points(cells, return_ids=True)
            cells = cells[ids]
            probs = probs[ids]

        # Turn off axes/grid
        scene.plotter.axes = 0  # 0 disables all axes and grids

        # vedo and brainrender have different coordinate systems, so align cells to brainrender HY
        cells[:, 2] *= -1

        points_actor_v = vPoints(cells, r=10)

        points_actor_v.cmap("viridis", probs, vmin=MIN_P)

        sc2d = ScalarBar(
            points_actor_v,
            #title="Membership probability",
            c="k",
            size=[0, 1],
            title_yoffset=20,
            font_size=20,
            #label_size=1.8,
            use_alpha = False,
            nlabels=5,
            label_format = "%.2f"
        )
        # colorscale formatting
        # pos
        sc2d.SetPosition(0.8, 0.05)   # center-bottom
        sc2d.SetWidth(0.13)
        sc2d.SetHeight(0.88)
        sc2d.SetTextPositionToSucceedScalarBar()   
        # labels
        label_prop = sc2d.GetLabelTextProperty()
        label_prop.SetFontSize(80)
        label_prop.SetBold(True)
        # title
        # Access the title text property # title behavior is super unexpected, so write it therem anually
        # title_prop = sc2d.GetTitleTextProperty()
        # title_prop.SetFontSize(40)
        # title_prop.SetBold(True)



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

        camera_settings = {
            "pos": pos,     # camera position in world coordinates
            "focalPoint": hy_com,     # point camera looks at
            "viewup": (0, -1, 0),         # "up" direction for the camera: -y ccf
        }

        outfile_name = f"{'_'.join(cluster_descr.split(" "))}_3D_cells_{camera_pos}_{MIN_P}"
        if RESTRICT_P_TO_HY:
            outfile_name += "_hyP"
        show(points_actor_v,sc2d, size=(2000, 2000), camera=camera_settings,
        screenshot=f"{str(OUTDIR)}/{outfile_name}.svg"
        )

        del scene

        # there may be a way to export this as an object, then to import into illustrator/blender as to get an svg