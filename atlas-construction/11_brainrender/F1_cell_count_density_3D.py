from pathlib import Path
import os 
import vedo
import numpy as np
import brainrender
from brainrender import Scene, settings
from brainrender.actors import Points,PointsDensity
import anndata as ad 
import pandas as pd 
from vedo import build_lut
from dotenv import dotenv_values

# settings
settings.ROOT_ALPHA = 0.05
settings.SHOW_AXES = False
vedo.settings.default_backend= 'vtk'
vedo.settings.default_font = "DejavuSansMono"

MAX_Y = 4000 

CONFIG = dotenv_values("../../config.env") 
DATADIR = Path(CONFIG['DATADIR'])

SPATIALDIR = DATADIR / "9_objects" / "spatial"
BRAINRENDER_DIR = DATADIR / "11_brainrender" 
OUTDIR = BRAINRENDER_DIR / "3D_plots"
os.makedirs(OUTDIR, exist_ok=True)

def get_camera(view="sagittal"):
    center = (6600, 4000, -5700)
    if view == "coronal":
        return dict(pos=(-28000, 4000, -5700), focalPoint=center, viewup=(0, -1, 0))
    elif view == "sagittal":
        return dict(pos=(6600, 4000, 33000), focalPoint=center, viewup=(0, -1, 0))
    elif view == "axial":
        return dict(pos=(6600, -45000, -5700), focalPoint=center, viewup=(0, 0, -1))
    if view == "frontleft":
        return dict(pos=(-17000, 4000, 25000), focalPoint=center, viewup=(0, -1, 0))


adata = ad.read_h5ad(SPATIALDIR / "hypomap3d_spatial_cell_level.h5ad", backed="r")

point_array = adata.obsm["coords_ccf_mm_wrap"]
point_array = np.array(point_array, dtype = "float32")
point_array = point_array * 1000 

# generate scene with hypo
scene = Scene(atlas_name='allen_mouse_25um', title=False,inset=False, root=True)

hy = scene.add_brain_region('HY', alpha=0.1, hemisphere="left"
    )

cut_plane = scene.atlas.get_plane(pos=(0, MAX_Y, 0), norm=(0, 1, 0))

scene.slice(cut_plane, actors=[hy])

dim_scale_factor = 1
actor_density = PointsDensity(point_array, 
    dims = (int(528 / dim_scale_factor), int(320 / dim_scale_factor), int(456 / dim_scale_factor)), #number of voxels in x, y and z of the output Volume. default = 40,40,40
    radius = 12.5 * dim_scale_factor, # our voxel size
    #radius = 100
    #compute_gradient=True #default false
    colors="viridis"
    )

# Build a custom Look-Up-Table of colors:
#     value, color, alpha
lut = build_lut(
    [
      #(-2, 'pink'      ),  # up to -2 is pink
      (0.9, "white", 0),
      #(1, 'white',  0.2    ),  # up to 0.9 is white with alpha = 0
      (1.9, 'blue', 0.2),  # up to 1.9 is blue with alpha=0.2
      (int(3 * dim_scale_factor), 'green', 0.7 ),
      (int(10 * dim_scale_factor), 'yellow', 0.7 ),
      #(int(10 * dim_scale_factor), 'green', 0.2 ),
      (int(25 * dim_scale_factor), 'white', 0.9 ), #25 is the max value we have

    ],
    interpolate=True,
)

radius = actor_density.metadata["radius"]

actor_density.cmap(lut)
actor_density.add_scalarbar(title=r"Cells per voxel", horizontal=False, nlabels=5) 

scalarbar_actor = actor_density.scalarbar

scene.add(actor_density,)


# Extract brainrender actors
brainrender_actors = scene.actors  # This is a list of vedo actors

all_actors = brainrender_actors + [scalarbar_actor]

#colorscale formatting
#vsettings.default_font = "DejavuSansMono"
#pos
scalarbar_actor.SetPosition(0.9, 0.3)   # center-bottom
scalarbar_actor.SetWidth(0.04)
scalarbar_actor.SetHeight(0.4)
scalarbar_actor.SetTextPositionToSucceedScalarBar()   
# labels
label_prop = scalarbar_actor.GetLabelTextProperty()
label_prop.SetFontSize(24)
label_prop.SetBold(True)
# title
title_prop = scalarbar_actor.GetTitleTextProperty()
title_prop.SetFontSize(24) 
title_prop.SetBold(True)

for camera_pos in ["coronal", "sagittal", "axial", "frontleft"]:
    camera = get_camera(camera_pos)  # or "axial", "sagittal"
    plot = vedo.Plotter()
    plot.show(*all_actors, axes=0, camera=camera, size=(1200, 800)
        )

    screenshot_path = OUTDIR / f"cells_per_voxel_3D_{camera_pos}.svg"
    plot.screenshot(screenshot_path)