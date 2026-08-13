"""
2_run_solo.py

This script runs the SOLO doublet detection model from scvi-tools on a subset of cells from a single-cell RNA-seq dataset,
using a pretrained scVI model. It is designed to be used within a SLURM pipeline and mode-specific filtering.

Usage:
    python 2_run_solo.py --version <VERSION> --mode <control|non_control> 

Arguments:
    --mode (str): 
        Which data subset to process. Must be one of:
        - "control": Uses `obs['IsControl'] == True`
        - "non_control": Uses `obs['IsControl'] == False`

Input Files:
    Located at: <DATA_ROOT_DIR>/3_solo/all_solo/datasets/
    - raw.h5ad: Full AnnData object with cell metadata
    - embed_only.h5ad: AnnData object containing reduced dimensions and graphs

Output:
    A CSV file of predicted doublets:
    <DATA_ROOT_DIR>/3_solo/all_solo/metrics/solo_doublets_<MODE>.csv

    - The file contains one column: `IsDoublet` for each cell in the relevant subset.

Notes:
    - Cells from technologies not requiring doublet detection (i.e. smart-seq) are excluded
    - Batches with fewer than 1024 cells are skipped (to avoid unstable SOLO training)
    - SOLO models are created from pretrained SCVI models and trained with GPU acceleration

Example:
    Run default control analysis:
        python 2_run_solo.py  --mode control

"""
# imports
import argparse
from pathlib import Path
import warnings

import scanpy as sc
import scvi
from torch import set_float32_matmul_precision
from dotenv import dotenv_values
import torch 

# seed for SCVI
scvi.settings.seed = 0

# Setup
set_float32_matmul_precision('high')
warnings.filterwarnings("ignore") #ignore view-related warnings from anndata

# Parse CLI arguments
parser = argparse.ArgumentParser()
parser.add_argument("--mode", type=str, choices=["control", "non_control"], required=True, help="Which part to process")
args = parser.parse_args()
MODE = args.mode

# Import config and construct/define variables
CONFIG = dotenv_values("../../config.env") # assuming kernel located at the root of the directory
DATA_ROOT_DIR = Path(CONFIG['DATADIR'])
MODELDIR = DATA_ROOT_DIR / "3_solo" / "all_solo"
BATCH_KEYS = ["Sample"]

# Specify technologies to exclude (e.g. SmartSeq datasets don't need doublet detection)
EXCLUDED_TECHNOLOGIES_REGEXES = [r"smart[-]?seq", ]
EXCLUDED_TECHNOLOGIES_PATTERN = "|".join(EXCLUDED_TECHNOLOGIES_REGEXES)

# helpers
def run_SOLO(vae, batch=None, **kwargs):
    if batch:
        model = scvi.external.SOLO.from_scvi_model(vae, restrict_to_batch=batch)
    else:
        model = scvi.external.SOLO.from_scvi_model(vae)
    model.train(**kwargs)
    return model

def get_batch_masks(adata, batch_keys):
    unique_combinations = adata.obs[batch_keys].drop_duplicates()
    for _, row in unique_combinations.iterrows():
        yield (adata.obs[batch_keys] == row.values).all(axis=1), "_".join(row.astype(str).values)


# Load data
print("Reading adata...")
adata = sc.read_h5ad(Path(MODELDIR, 'datasets', 'raw.h5ad'))
embed_only = sc.read_h5ad(Path(MODELDIR, 'datasets', 'embed_only.h5ad'), backed='r')

adata.obsm = embed_only[adata.obs_names].obsm
adata.uns = {k: v for k, v in embed_only.uns.items() if k in {'tsne', 'umap', 'neighbors'}}
adata.obsp = embed_only.obsp
del embed_only

# init IsDoublet
adata.obs['IsDoublet'] = False
# restrict to valid technologies
valid_mask = ~(adata.obs["Technology"].astype(str).str.lower().str.contains(EXCLUDED_TECHNOLOGIES_PATTERN, na=False))
valid_subset = adata[valid_mask]
print(f"{len(valid_subset)} valid cells")

if MODE == "control":
    adata_subset = valid_subset[valid_subset.obs['IsControl']]
elif MODE == "non_control":
    adata_subset = valid_subset[~valid_subset.obs['IsControl']]

for mask, name in get_batch_masks(adata_subset, BATCH_KEYS):
    n_cells = sum(mask)
    print(f"Current batch: {name} ({n_cells} cells)")
    batch_subset = adata_subset[mask].copy()
    if n_cells >= 1024:
        scvi_model_dir = Path(MODELDIR, 'models', 'scvi').as_posix()
        scvi.model.SCVI.prepare_query_anndata(batch_subset, scvi_model_dir)
        scvi_model = scvi.model.SCVI.load(scvi_model_dir, adata=batch_subset)
        solo_model = run_SOLO(scvi_model,
                                batch=name,
                                accelerator="gpu" if torch.cuda.is_available() else "cpu",
                                max_epochs=400,
                                batch_size=1024,
                                datasplitter_kwargs={'drop_last': True})
        adata.obs.loc[batch_subset.obs_names, 'IsDoublet'] = solo_model.predict(soft=False) == 'doublet'
    else:
        print('\t - Too small to run SOLO')
        adata.obs.loc[batch_subset.obs_names, 'IsDoublet'] = False

# Save only the relevant subset based on mode
if MODE == "control":
    output_subset = adata[adata.obs['IsControl']]
else:
    output_subset = adata[~adata.obs['IsControl']]

output_path = Path(MODELDIR, 'metrics')
output_path.mkdir(exist_ok=True, parents=True)

output_file = output_path / f'solo_doublets_{MODE}.csv'
output_subset.obs[['IsDoublet']].to_csv(output_file)
print(f"Saved predictions to {output_file}")

