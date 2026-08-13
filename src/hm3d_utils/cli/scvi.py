"""
This script sets up and runs an integration model for single-cell RNA-seq data analysis using scVI. 
It includes steps to load and preprocess the dataset, train a scVI model, and compute dimensionality 
reduction embeddings such as UMAP and t-SNE. 

Overview:
- `main`: Main function to execute the integration pipeline.
- `initialize_model_directory`: Creates a directory structure for storing model results, parameters, and datasets.
- `get_model_embedding`: Retrieves the latent representation of cells using a specified model function.
- `get_scvi_embedding`: Trains or loads a pre-trained scVI model to obtain the latent embedding.
- `plot_model_history`: Plots the training history of the scVI model, including ELBO and reconstruction loss.
  
Dependencies:
- argparse: For parsing command-line arguments.
- json: For loading model parameters from a JSON file.
- logging: For logging informational and debugging messages.
- os: For file and directory manipulation.
- shutil: For copying files.
- warnings: To ignore deprecation and user warnings.
- numpy, pandas: For numerical and data manipulation.
- scanpy: For processing and analyzing single-cell RNA-seq data.
- scvi: For training and using scVI models.
- torch: For setting floating point precision in matrix multiplications.
- matplotlib: For plotting the model's training history.
- rapids_singlecell: (Optional) For GPU-accelerated computation.

Functions:
- main(args: argparse.Namespace) -> None: Main pipeline execution function, orchestrating data loading, model training, 
  and embedding calculations.
- initialize_model_directory(parameter_file: str) -> None: Creates directory structure and copies input files for 
  the integration model.
- get_model_embedding(adata: sc.AnnData, parameters: dict, model_func) -> np.array: Obtains the latent embedding 
  using a provided model function.
- get_scvi_embedding(adata: sc.AnnData, parameters: dict) -> np.array: Loads or trains a scVI model to get the 
  latent representation of the input data.
- plot_model_history(history) -> matplotlib.figure.Figure: Generates and returns a plot of the scVI model's training history.

Example usage:
    python script.py integration_params.json --verbose
"""

import argparse
import json
import logging
import os
from shutil import copy2, SameFileError
import warnings
from pathlib import Path

logger = logging.getLogger(__name__)
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from torch import set_float32_matmul_precision
from dotenv import dotenv_values 

scvi.settings.seed=0
set_float32_matmul_precision('medium')

global PARAMETERS
global OUTDIR
EMBEDDING_KEY = 'X_embed'

def main(args: argparse.Namespace) -> None:
    """
    Main function to execute the single-cell RNA-seq integration pipeline using scVI. 
    This function initializes model directories, loads datasets, subsets genes, and calculates embeddings, 
    including UMAP, t-SNE, and Leiden clustering.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments containing the path to the integration parameter JSON file 
        and optional verbosity level.

    Returns
    -------
    None
    """
    
    global PARAMETERS
    global OUTDIR 
    
    with open(args.parameter_file) as parameter_file:
        PARAMETERS = json.load(parameter_file)

    #prepend data root dir to paths
    CONFIG = dotenv_values("../../config.env")
    root = Path(CONFIG["DATADIR"])
    for k in ["outdir", "hvgs", "dataset"]:
        PARAMETERS[k] = str(root / PARAMETERS[k])
    
    OUTDIR = os.path.join(PARAMETERS['outdir'], PARAMETERS['model_name'])

    logging.info(f'Initializing {OUTDIR}')
    
    initialize_model_directory(args.parameter_file)
    
    adata = sc.read(PARAMETERS['dataset'])
    logger.info(f'Successfully loaded raw dataset.')
    
    hvgs = pd.read_csv(PARAMETERS['hvgs'], index_col=0, header=None)
    logger.info(f'Subsetting to highly variable genes ({len(hvgs)}).')
    adata = adata[:, hvgs.index].copy()

    logger.info(adata)

    logger.info('Getting model embedding.')
    adata.obsm[EMBEDDING_KEY] = get_model_embedding(adata, PARAMETERS, get_scvi_embedding)

    del adata.X

    logger.info('Calculating neighbors.')
    neighbors_parameters = PARAMETERS.get('neighbors', {})
    sc_neighbors_kwargs = {"use_rep": EMBEDDING_KEY} | neighbors_parameters.get('kwargs', {})
    rsc_neighbors_kwargs = {"use_rep": EMBEDDING_KEY} | neighbors_parameters.get('kwargs', {})
    if neighbors_parameters.get('gpu', False):
        try:
            import rapids_singlecell as rsc
            rsc.pp.neighbors(adata, **rsc_neighbors_kwargs)
        except ModuleNotFoundError as e:
            logger.warning(f'RAPIDS gpu accelration unavailable, falling back to cpu. \n\t!> {e}')
            sc.pp.neighbors(adata, **sc_neighbors_kwargs)   
    else:
        sc.pp.neighbors(adata, **sc_neighbors_kwargs)
        
    logger.info('Calculating UMAP.')
    umap_parameters = PARAMETERS.get('umap', {})
    sc_umap_kwargs = {} | umap_parameters.get('kwargs', {})
    rsc_umap_kwargs = {} | umap_parameters.get('kwargs', {})
    if umap_parameters.get('gpu', False):
        try:
            import rapids_singlecell as rsc
            rsc.tl.umap(adata, **rsc_umap_kwargs)
        except ModuleNotFoundError as e:
            logger.warning(f'RAPIDS gpu accelration unavailable, falling back to cpu. \n\t!> {e}')
            sc.tl.umap(adata, **sc_umap_kwargs)
    else:
        sc.tl.umap(adata, **sc_umap_kwargs)
        
    logger.info('Calculating t-SNE.')
    tsne_parameters = PARAMETERS.get('tsne', {})
    sc_tsne_kwargs = {"perplexity": 10} | tsne_parameters.get('kwargs', {})
    rsc_tsne_kwargs = {"perplexity": 10} | tsne_parameters.get('kwargs', {})
    if tsne_parameters.get('gpu', False):
        try:
            import rapids_singlecell as rsc
            rsc.tl.tsne(adata, use_rep=EMBEDDING_KEY, **rsc_tsne_kwargs)
        except ModuleNotFoundError as e:
            logger.warning(f'RAPIDS gpu accelration unavailable, falling back to cpu. \n\t! {e} !')
            sc.tl.tsne(adata, use_rep=EMBEDDING_KEY, **sc_tsne_kwargs)
    else:
        sc.tl.tsne(adata, use_rep=EMBEDDING_KEY, **sc_tsne_kwargs)
        
#    logger.info('Calculating leiden.')
#    leiden_parameters = PARAMETERS.get('leiden', {})
#    sc_leiden_kwargs = {"flavor": "igraph", "n_iterations": 2} | leiden_parameters.get('kwargs', {})
#    rsc_leiden_kwargs = {} | leiden_parameters.get('kwargs', {})
#    if leiden_parameters.get('gpu', False):
#        try:
#            import rapids_singlecell as rsc
#            rsc.tl.leiden(adata, **rsc_leiden_kwargs)
#        except ModuleNotFoundError as e:
#            print(f'RAPIDS gpu accelration unavailable, falling back to cpu. \n\t! {e} !')
#            sc.tl.leiden(adata, **sc_leiden_kwargs)
#    else:
#        sc.tl.leiden(adata, **sc_leiden_kwargs)
    
    embedding_only_file = os.path.join(OUTDIR, 'datasets/embed_only.h5ad')
    adata.write_h5ad(embedding_only_file)
    logger.info(f'Saving embedding to {embedding_only_file}')

def initialize_model_directory(parameter_file: str) -> None:
    """
    Create directory for new model (WORKDIR), hard linking datasets and copying the parameter file.

    Parameters
    ----------
    parameter_file : str
        Filepath to the parameter file.

    Returns
    -------
    None
    """
    global PARAMETERS
    global OUTDIR
    
    os.makedirs(os.path.join(OUTDIR, 'datasets'), exist_ok=True)
    os.makedirs(os.path.join(OUTDIR, 'metrics'), exist_ok=True)
    os.makedirs(os.path.join(OUTDIR, 'models'), exist_ok=True)
    os.makedirs(os.path.join(OUTDIR, 'plots'), exist_ok=True)
    logger.info(f'Subdirectories created in {OUTDIR}')

    try:
        copy2(PARAMETERS['hvgs'], os.path.join(OUTDIR, 'metrics/hvgs.csv'))
        copy2(PARAMETERS['dataset'], os.path.join(OUTDIR, 'datasets/raw.h5ad'))
        logger.info(f'Copied datasets into {OUTDIR}')
    except SameFileError:
        logger.warn(f'Datasets already present in {OUTDIR}')
        pass
    
    try:
        copy2(parameter_file, os.path.join(OUTDIR,'parameters.json'))
    except SameFileError:
        logger.warn(f'Parameter file is already present in {OUTDIR}')
        pass

def get_model_embedding(adata: sc.AnnData, parameters: dict, model_func) -> np.array:
    """
    Obtains the latent embedding of cells using the specified model function.

    Parameters
    ----------
    adata : sc.AnnData
        Annotated data matrix from the `scanpy` library.
    
    parameters : dict
        Dictionary containing model parameters.
    
    model_func : function
        A function to compute the model embedding.

    Returns
    -------
    np.array
        Array of the latent embeddings.
    """
    return model_func(adata, parameters)

def get_scvi_embedding(adata: sc.AnnData, parameters: dict) -> np.array:  
    """
    Trains or loads a scVI model and obtains the latent representation of the input AnnData object.

    Parameters
    ----------
    adata : sc.AnnData
        Annotated data matrix from the `scanpy` library.
    
    parameters : dict
        Dictionary containing parameters for setting up and training the scVI model.

    Returns
    -------
    np.array
        The latent representation of the cells in the AnnData object.
    """
    global OUTDIR
    
    model_dir = os.path.join(OUTDIR, 'models/scvi')
    model_parameters = parameters['model']
    logger.info(model_parameters)
    
    if os.path.exists(model_dir):
        model = scvi.model.SCVI.load(model_dir, adata=adata)
    else:
        scvi.model.SCVI.setup_anndata(adata, **model_parameters.get('setup_anndata_kwargs', {}))
        model = scvi.model.SCVI(adata, **model_parameters.get('model_kwargs', {}))
        model.train(**model_parameters.get('train_kwargs', {}))
        model.save(model_dir)
    
    fig = plot_model_history(model.history)
    fig.savefig(os.path.join(OUTDIR, 'plots/model_history.png'))

    return model.get_latent_representation()

def plot_model_history(history):
    """
    Plots the training history of the scVI model, including ELBO and reconstruction loss over epochs.

    Parameters
    ----------
    history : dict
        A dictionary containing the training history of the model, with keys such as 'elbo_train', 
        'elbo_validation', 'reconstruction_loss_train', 'reconstruction_loss_validation', and 'kl_weight'.

    Returns
    -------
    matplotlib.figure.Figure
        A figure object containing the plotted training history.
    """
    import matplotlib.pyplot as plt
    
    fig = plt.figure(figsize=(10, 20))
    gs = fig.add_gridspec(3, hspace=0)
    axs = gs.subplots(sharex=True)

    axs[0].plot(history['elbo_train'], label='Train', color='black')
    axs[0].plot(history['elbo_validation'], label='Validation', color='red')
    axs[0].set_ylabel('ELBO')
    axs[0].legend()

    axs[1].plot(history['reconstruction_loss_train'], label='Train', color='black')
    axs[1].plot(history['reconstruction_loss_validation'], label='Validation', color='red')
    axs[1].set_ylabel('Reconstruction Loss')
    axs[1].legend()

    axs[2].plot(history['kl_weight'], label='Training KL Weight', color='purple')
    axs[2].set_ylabel('KL Weight')
    axs[2].set_xlabel('Epoch')
    axs[2].xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    axs[2].legend()

    for ax in axs:
        ax.label_outer()
        
    plt.close()

    return fig

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Setup raw data for integration model setup')

    parser.add_argument('parameter_file', metavar='[ Integration Parameters (JSON) ]', help='Integration Parameters')
    parser.add_argument('--verbose', '-v', action='count', default=0, help='Increase output verbosity (e.g., -vv for more)')
    
    args = parser.parse_args()

    log_level = {0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}
    logging.basicConfig(level=log_level.get(args.verbose, logging.DEBUG), 
                        format='%(asctime)s - %(levelname)s - %(message)s', 
                        datefmt='%Y-%m-%d %H:%M:%S')

    main(args)
