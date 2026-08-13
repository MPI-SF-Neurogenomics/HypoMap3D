from pathlib import Path
import argparse
import os

import pandas as pd
import scanpy as sc
from dotenv import dotenv_values

from hm3d_utils.statistics import find_markers
from hm3d_utils.utils import get_clusters_in_region


############## ############## ##############
############## CLI arguments
############## ############## ##############

parser = argparse.ArgumentParser(
    description="Calculate region-specific marker genes."
)

parser.add_argument(
    "--include-nn",
    action="store_true",
    help="Include NN cell types."
)

parser.add_argument(
    "--analysis-level",
    type=int,
    choices=[1, 2, 3, 4, 5],
    default=3,
    help="Annotation level (L1-L5). Default: 3"
)

parser.add_argument(
    "--acro",
    type=str,
    required=True,
    help="Region acronym (e.g. DMH, PVH, ARH)"
)

args = parser.parse_args()

INCLUDE_NN = args.include_nn
ANALYSIS_LEVEL = args.analysis_level
ACRO = args.acro


############## ############## ############## 
############## Environment
############## ############## ############## 

CONFIG = dotenv_values("../../config.env") 
DATADIR = Path(CONFIG['DATADIR'])

# adata
ADATA_FP = DATADIR / "6_NC_projection" / "all_integrated" / "HypoMap3D_full.h5ad"

# regions
REGIONS_ALL_FP = DATADIR / "7_spatial" / "analyses" / "regionwise_rm_all_wo_root_or_fibers.csv" # allowing all regions and HY acronym

# out
OUTDIR = DATADIR / "8_annotation" / "region_specific_markers"
OUTDIR.mkdir(parents=True, exist_ok=True) 

############## ############## ############## 
############## OUTFILE handling
############## ############## ############## 

if not INCLUDE_NN:
    OUTFILE_MARKERS = OUTDIR / f"L{ANALYSIS_LEVEL}_markers_{ACRO}_specific.csv"
    OUTFILE_CLUSTERS = OUTDIR / f"L{ANALYSIS_LEVEL}_{ACRO}_clusters.csv"
else:
    OUTFILE_MARKERS = OUTDIR / f"L{ANALYSIS_LEVEL}_markers_{ACRO}_specific_w_NN.csv"
    OUTFILE_CLUSTERS = OUTDIR / f"L{ANALYSIS_LEVEL}_{ACRO}_clusters_w_NN.csv"

if os.path.exists(OUTFILE_MARKERS):
    os.remove(OUTFILE_MARKERS)

############## ############## ############## 
############## Data loading and preprocessing 
############## ############## ############## 

# adata
adata = sc.read_h5ad(ADATA_FP, backed="r")

# subset adata 
load_mask = adata.obs.IsControl
if not INCLUDE_NN:
    load_mask = load_mask & (~adata.obs.L1.str.contains("NN"))

adata = adata[load_mask,:].to_memory() # subset to control cells
sc.pp.normalize_total(adata, target_sum = 10_000) # run markers function operates on non-log CP10K


# region predictions 
regions = pd.read_csv(REGIONS_ALL_FP)
regions = regions.rename(columns={"Unnamed: 0": "cluster"})

# filter to cell types allowed in the analysis 
load_mask = regions.cluster.str[2] == str(ANALYSIS_LEVEL)
if not INCLUDE_NN:
    load_mask = load_mask & (~regions.cluster.str.contains("NN"))

regions = regions.loc[load_mask, :]
regions = pd.melt(regions, id_vars="cluster", value_name="relative_mass", var_name = "acronym")
regions = regions.sort_values(by=["cluster", "relative_mass"], ascending=False)

############## ############## ############## 
############## Marker calculation 
############## ############## ############## 

cluster_df = get_clusters_in_region(regions, acronym=ACRO, top_k=5, min_rel_mass=0.1, cum_mass_max=0.6, cum_mass_mode="excluding")
cluster_df.to_csv(OUTFILE_CLUSTERS, index=False)
clusters = cluster_df.cluster.unique()
n_clusters = len(clusters)
write_header = True

if n_clusters < 2:
    print(f"Less than 2 clusters found. Skipping markers calculation.")
else:
    print(f"{n_clusters} clusters found. Calculating markers...")

    for i, cluster in enumerate(clusters):

        print(f"Processing {i+1}/{n_clusters}")

        # get masks
        others = [c for c in clusters if c != cluster]

        self_mask = adata.obs[f"L{ANALYSIS_LEVEL}"] == cluster
        others_mask = adata.obs[f"L{ANALYSIS_LEVEL}"].isin(others)

        try:
            markers_df = find_markers(
                adata_A=adata[self_mask],
                adata_B=adata[others_mask],
                strata_key="Sample",
                max_stratum_size=10_000,
                min_stratum_size=10,
                min_pct_expression=0.1,
                min_specificity=0.1,
            ).reset_index(names = "gene_name")

            # filter
            markers_df = markers_df.loc[(markers_df.direction == "A") & (
                markers_df.log2FC > 0) & (markers_df.padj < 0.05
            ), [c for c in markers_df.columns if c != "direction"]]

            # add metadata
            markers_df["region"] = ACRO
            markers_df["cluster"] = cluster

            # flush to csv
            markers_df.to_csv(
                OUTFILE_MARKERS,
                mode="a",
                header=write_header,
                index=False
            )

            write_header = False

            print(f"{cluster}: {len(markers_df)} markers found.")
            
        except Exception as e:
            print(f"    Warning: run_markers failed for cluster '{cluster}' in {ACRO}: {e}", flush=True)
            continue
