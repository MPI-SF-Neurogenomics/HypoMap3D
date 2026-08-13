"""
This script serves to precompute the coordinate/voxel space for 
faster access when plotting voxel-level data.

The script assumes your that you are running it with its directory being 
your working directory.

Then, you can run 

singularity exec path/to/uv/image uv run 0_precompute_coord_space.py
"""

import os 
import zarr
from pathlib import Path

import pandas as pd
import numpy as np 
from brainrender import Scene
from dotenv import dotenv_values


CONFIG = dotenv_values("../../config.env") 
DATADIR = Path(CONFIG['DATADIR'])

SPATIAL_DIR = DATADIR / "9_objects" / "spatial"
OUTDIR = DATADIR / "11_brainrender" 
os.makedirs(OUTDIR, exist_ok=True)
OUTFILE = OUTDIR / "coord_space.csv"

ZARR_FP = SPATIAL_DIR / "hypomap3d_spatial_voxel_level.zarr"
REGION_META_FP = DATADIR / "7_spatial" / "abc_region_metadata.csv"

MAX_Y = 4000 # used to remove the SFO


# load spatial data
voxel_store = zarr.open(ZARR_FP)


# get flat voxel grid
def get_bbox_ccf_coords_flat(zarrfile):

    min_coords = zarrfile["bbox"]["ccf_mm"]["min"] # min_x, min_y, min_z, shape (3,)
    max_coords = zarrfile["bbox"]["ccf_mm"]["max"]

    xx = np.arange(start=min_coords[0], stop=max_coords[0]+0.025, step=0.025) 
    yy = np.arange(start=min_coords[1], stop=max_coords[1]+0.025, step=0.025) 
    zz = np.arange(start=min_coords[2], stop=max_coords[2]+0.025, step=0.025) 

    # Create 3D meshgrid from coordinate arrays
    X, Y, Z = np.meshgrid(xx, yy, zz, indexing='ij')

    # Stack into a single array of shape (nx, ny, nz, 3)
    grid = np.stack([X, Y, Z], axis=-1)

    coords = grid.reshape(-1, 3)

    return coords

df = pd.DataFrame(data=get_bbox_ccf_coords_flat(voxel_store), columns=["x_ccf", "y_ccf", "z_ccf"])
df.x_ccf *= 1000 # transform from mm to um
df.y_ccf *= 1000
df.z_ccf *= 1000

# add allen data
region = voxel_store["region"]
# Flatten voxel-wise Allen IDs
df["allen_id"] = region["region_labels"][:].ravel()
# Create ID -> acronym mapping
region_meta = pd.read_csv(REGION_META_FP)
id_to_acronym = region_meta.set_index("ID").to_dict()["acronym"]
# Map IDs to acronyms
df["allen_acronym"] = df["allen_id"].map(id_to_acronym)

# Create scene with the HY 
scene = Scene(root=False)
hy = scene.add_brain_region('HY', alpha=0.1, color="red", hemisphere="right"
    )
# remove the SFO
cut_plane = scene.atlas.get_plane(pos=(0, MAX_Y, 0), norm=(0, 1, 0))

scene.slice(cut_plane, actors=[hy])


# check what coordinates are within the HY hull
coords = df[['x_ccf', 'y_ccf', 'z_ccf']].values
ids = hy.mesh.inside_points(coords, return_ids=True)

# insert indicator 
df["is_inside_HY"] = False
df.loc[ids, "is_inside_HY"] = True

# save
df.to_csv(OUTFILE,index=False)
print(f"File created at {OUTFILE}")
