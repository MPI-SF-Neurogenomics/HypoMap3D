import os
import gzip
import shutil
import multiprocessing
from pathlib import Path
import argparse
import logging

def setup_logger(log_dir: Path):
    """
    Set up logging to both console and a file located at log_dir/compression.log.
    Creates the log directory if it doesn't exist.
    """
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / 'compression.log'

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

def gzip_file(fastq_path: Path):
    """
    Compress a .fastq file using gzip and validate the result.
    Skips compression if the .gz version already exists.
    Logs all actions and errors.
    """
    gz_path = fastq_path.with_suffix(fastq_path.suffix + '.gz')

    # Skip if gzipped version already exists
    if gz_path.exists():
        logging.info(f"Skipping (already gzipped): {fastq_path.name}")
        return

    try:
        # Compress the file
        logging.info(f"Compressing: {fastq_path.name} -> {gz_path.name}")
        with open(fastq_path, 'rb') as f_in, gzip.open(gz_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

        # Validate the compressed file by trying to read 1 byte
        with gzip.open(gz_path, 'rb') as f_check:
            f_check.read(1)
        logging.info(f"Validated: {gz_path.name} is readable and appears intact.")

    except Exception as e:
        logging.error(f"Failed during compression/validation of {fastq_path.name}: {e}")

def main(directory: str):
    """
    Main function that:
    - Validates the input directory
    - Sets up logging
    - Finds all .fastq files without a .gz version
    - Compresses them using multiprocessing
    """
    dir_path = Path(directory).resolve()
    if not dir_path.is_dir():
        print(f"Error: {dir_path} is not a valid directory.")
        return

    # Create log dir
    log_dir = dir_path / 'logs'
    setup_logger(log_dir)

    logging.info(f"Starting compression in directory: {dir_path}")

    # Gather all .fastq files that are not yet gzipped
    fastq_files = [
        f for f in dir_path.iterdir()
        if f.is_file() and f.suffix == '.fastq' and not (f.with_suffix(f.suffix + '.gz')).exists()
    ]

    if not fastq_files:
        logging.info("No .fastq files to compress.")
        return

    logging.info(f"Found {len(fastq_files)} .fastq files to compress.")

    # Compress files in parallel
    with multiprocessing.Pool() as pool:
        pool.map(gzip_file, fastq_files)

    logging.info("Compression complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Gzip all .fastq files in a directory if not already gzipped. "
                    "Validates each compressed file. Logs actions to /logs."
    )
    parser.add_argument("directory", help="Path to the directory containing .fastq files.")
    args = parser.parse_args()

    main(args.directory)
