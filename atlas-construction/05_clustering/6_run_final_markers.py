from pathlib import Path

import scanpy as sc
from dotenv import dotenv_values

from hm3d_utils.statistics import get_final_markers
from hm3d_utils.utils import ClusterTree


def main():

    CONFIG = dotenv_values("../../config.env") 
    DATADIR = Path(CONFIG['DATADIR'])
    OUTDIR = DATADIR / "5_clustering" 

    ADATA_FP = OUTDIR / "all_control" / "datasets" / 'control_clustered.h5ad'
    OUTFILE_PATH = OUTDIR / 'HypoMap3D_marker_genes.csv'

    print("Loading AnnData...")
    adata = sc.read_h5ad(ADATA_FP)

    print("Normalizing gene expression...") 
    sc.pp.normalize_total(adata, target_sum=10_000)
    adata.X.eliminate_zeros()

    print("Creating trees...")
    nn_mask = adata.obs["L1"].str.contains("NN")

    n_cluster_tree = ClusterTree.from_labelmat(
                                adata.obs.loc[~nn_mask, [f"L{i}" for i in range(1,6)]].values, # L1-5
                                adata[~nn_mask].obsm["X_embed"], 
                                name="Neuronal") 

    nn_cluster_tree = ClusterTree.from_labelmat(
                                adata.obs.loc[nn_mask, [f"L{i}" for i in range(1,5)]].values, # L1-4
                                adata[nn_mask].obsm["X_embed"], 
                                name="Non-neuronal")

    print("Creating lookup tables...")
    n_lookup = n_cluster_tree.generate_lookup_table()
    n_lookup.pop("Neuronal") 
    nn_lookup = nn_cluster_tree.generate_lookup_table()
    nn_lookup.pop("Non-neuronal") 

    get_final_markers(
        adata=adata, 
        n_lookup=n_lookup, 
        nn_lookup=nn_lookup, 
        output_file=OUTFILE_PATH
    )

if __name__ == "__main__":
    main()
