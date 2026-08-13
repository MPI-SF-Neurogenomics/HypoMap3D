import argparse
import scanpy as sc
import pandas as pd

def main(args: argparse.Namespace) -> None:
    # Load the 10X Genomics data from the HDF5 file
    adata = sc.read_10x_h5(args.h5)
    print(f'Loaded data from {args.h5}')

    # If barcodes file is provided, filter the AnnData object
    if args.barcodes:
        selected_barcodes = pd.read_csv(args.barcodes, header=None)[0].tolist()
        adata = adata[adata.obs_names.isin(selected_barcodes)]
        print(f'Filtered to {len(adata)} cells using barcodes from {args.barcodes}')
    
    # Write the AnnData object to the specified output file
    adata.write_h5ad(args.outfile)
    print(f'Saved AnnData to {args.outfile}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert 10X .h5 data to .h5ad format with optional barcode filtering')

    parser.add_argument('h5', metavar='H5_FILE', type=str,
                        help='Path to the 10X Genomics .h5 file (filtered or raw)')
    parser.add_argument('-b', '--barcodes', type=str, default=None,
                        help='Path to a CSV file containing barcodes to filter the data')
    parser.add_argument('-O', '--outfile', type=str, default='./output.h5ad',
                        help='Path to the output .h5ad file')

    args = parser.parse_args()
    main(args)
