import anndata as ad 
import pandas as pd 
import re 

def standardize_obs_names(
    adata: ad.AnnData,
    dataset_name: str,
    sample_name: str,
    add_info_to_obs: bool = True,
):
    """
    Standardize cell barcodes / observation names across datasets.
    """

    # Convert obs_names to string just in case
    obs_names = adata.obs_names.astype(str)

    # 1) Remove trailing "-1" Cell Ranger artifact
    #    e.g. GTATTGGAGCCGATCC-1 -> GTATTGGAGCCGATCC
    obs_names = obs_names.str.replace(r"-1$", "", regex=True)

    # 2) Smart-seq special case:
    #    use SRR accession as cell ID
    #    e.g. '../data/v2-04-SRR26100756' -> 'SRR26100756'.
    #    If the first obs_name contains an SRR pattern,
    #    extract SRR IDs for all obs_names
    first_name = obs_names[0]
    if pd.notna(first_name) and re.search(r"SRR\d+", first_name):
        obs_names = obs_names.str.extract(r"(SRR\d+)", expand=False)

    # Update obs_names
    obs_names = dataset_name + "_" + sample_name + "_" + obs_names
    adata.obs_names = obs_names

    # 3) Add info to obs
    if add_info_to_obs:
        adata.obs["CellID"] = adata.obs_names
        adata.obs["Dataset"] = dataset_name
        adata.obs["Sample"] = sample_name

