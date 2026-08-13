import argparse
import logging
# import warnings
import sys

from pathlib import Path

import anndata as ad
import pandas as pd

# warnings.filterwarnings('ignore', category=DeprecationWarning)
# warnings.filterwarnings("ignore", category=UserWarning)

def main(args):
    import scvi
    
    from torch import set_float32_matmul_precision
    from torch.cuda import is_available
    
    scvi.settings.seed = args.seed
    set_float32_matmul_precision('medium')
   
    # Setup logging (to STDERR)
    verbose = args.verbose
    logger = setup_logger(verbose)
    
    model_path = Path(args.model)
    query_path = Path(args.query)
    outfile = Path(args.outfile)
    
    # Pre-flight check to ensure final write directory exists
    if not outfile.parent.exists():
        raise FileNotFoundError(f"Directory `{outfile.parent}` does not exist.")

    # Reformat query AnnData to match reference AnnData structure
    logger.info("Loading and formatting Query AnnData.")
    qry_adata = ad.read_h5ad(query_path)
    scvi.model.SCVI.prepare_query_anndata(qry_adata, model_path.as_posix(), inplace=True)
    
    logger.info(qry_adata)
    
    # Load the model with scArches surgery for new batches
    logger.info("Loading reference model with Query modifications.")
    qry_vae = scvi.model.SCVI.load_query_data(
        qry_adata,
        model_path.as_posix()
    )
    
    # Re-train with frozen hidden weights
    logger.info("Beginning training with Query AnnData.")
    qry_vae.train(
        accelerator="gpu" if is_available() else "cpu",
        max_epochs=args.epochs, 
        batch_size=args.batch_size,
        load_sparse_tensor=True,
        check_val_every_n_epoch=1,
        plan_kwargs=dict(
            weight_decay=0.0 # important for scArches
        ),
        datasplitter_kwargs=dict(
            drop_last=True
        ),
        early_stopping=True,
        early_stopping_patience=4,
        # early_stopping_warmup_epochs=8, # not supported in 1.2.2 :(
    )
    
    # Extract embedding from model
    logger.info("Extracting query embedding.")
    qry_X_embed = qry_vae.get_latent_representation()
    output = pd.DataFrame(
        qry_X_embed,
        columns=[f"latent_{i}" for i in range(qry_X_embed.shape[1])],
        index=qry_adata.obs.index
    )
    
    # Save to disk embeddings
    logger.info(f"Saving query embedding to disk - {outfile}")
    compress = outfile.suffixes[-2:] == ['.csv', '.gz']
    compression = "gzip" if compress else None
    output.to_csv(outfile, compression=compression)
    
def setup_logger(verbose=False):
    
    level = logging.INFO if verbose else logging.WARNING
    
    logger = logging.getLogger(__name__)
    logger.setLevel(level)

    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(level)

    formatter = logging.Formatter(
        "%(levelname)s: %(message)s"
    )
    handler.setFormatter(formatter)

    logger.handlers.clear()      # avoid duplicate handlers
    logger.addHandler(handler)
    logger.propagate = False
    
    return logger

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
        description='Projection of new query data onto a precomputed scVI model via scArches.'
    )
    
    # Positional Arguments
    parser.add_argument(
        'model',
        help='Path to directory containing saved model weights (`model.pt`)'
    )
    parser.add_argument(
        'query',
        help='Filepath to Query AnnData (.h5ad)'
    )
    parser.add_argument(
        "-o", "--outfile", 
        help="CSV file to which to write output embedding results.",
        required=True
    )
    
    # Optioanl Arguments
    parser.add_argument(
        "--epochs", 
        type=int, default=20,
        help="Number of training epochs."
    )
    parser.add_argument(
        "--batch_size", 
        type=int, default=128,
        help="Mini-batch size for training."
    )
    parser.add_argument(
        "--seed", 
        type=int, default=0,
        help="Random seed for reproducibility."
    )
    parser.add_argument(
        "-v", "--verbose", 
        action="store_true",
        help="Adding this flag will show INFO level logging information from process to STDERR."
    )
    
    args = parser.parse_args()
    
    main(args)