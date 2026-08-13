import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import scanpy as sc

def normalize_smartseq(adata: sc.AnnData, path_to_transcript_lengths:str) -> sc.AnnData:

    print("Normalizing gene lengths for Smart Seq data correction")

    gene_len = pd.read_csv(
    path_to_transcript_lengths,
    delimiter=",",
    header=None,
    index_col=0,
    )

    gene_len.index.name = 'gene_ids'
    gene_len = gene_len.rename(columns={1: 'canonical_transcript_length'})

    var_temp = adata.var.join(gene_len, on='gene_ids')

    adata.X = adata.X.toarray() / var_temp['canonical_transcript_length'].values * np.median(var_temp['canonical_transcript_length'].values)
    adata.X = csr_matrix(np.rint(adata.X))

    return adata