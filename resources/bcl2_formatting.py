import os
import glob
import re
import argparse
import pandas as pd

import multiprocessing
import logging
from logging.handlers import QueueHandler

from hm3d_utils.logging import (
    create_shared_logger_process_for_multiprocessing,
    get_shared_queue_and_logger_instances,
)

from hm3d_utils.formatting import (
    get_SRR_name_SRR_I_fastqgz, 
    process_srr_to_sample_symlink,
    process_bcl2_to_sample_symlink,  
    get_log_file_path_symlink,
    get_filename_stem,
)

global DATADIR
# some datasets don't have S info, so S info will be inferred later when processing these bcl2-like files.
BCL2_PATTERN = r'L[0-9]+[_\-.][RI][0-9]+[_\-.][0-9]+\.fastq\.gz' # some datasets don't have S info


def main_process_srr_to_sample(args: argparse.Namespace) -> None:
    """
    Main function for processing FASTQ files, mapping SRR IDs to sample names, and creating symbolic links
    for the FASTQ files. This function handles both standard FASTQ files and BCL2-related files separately,
    and uses multiprocessing to handle the workload efficiently. It also sets up logging for subprocesses.

    The function performs the following steps:
    1. Collects all `.fastq.gz` files from the specified directory.
    2. Separates the files into two categories: BCL2-related files and non-BCL2 files, based on pattern matching.
    3. Extracts SRR IDs from non-BCL2 files and the base filenames from BCL2 files.
    4. Sets up a shared logger to handle logging across multiple worker processes.
    5. Launches separate worker processes for processing non-BCL2 files and BCL2-related files.
    6. Waits for all worker processes to finish and ensures that the logger process terminates properly.

    Args:
        args (argparse.Namespace): Command-line arguments passed to the script, containing at least the `dir`
                                    argument, which specifies the directory containing the FASTQ files to be processed.

    Raises:
        Exception: If an error occurs during the execution of the worker processes or logger process.
    
    Notes:
        - This function uses the `multiprocessing` module to process files in parallel, improving performance.
        - It processes non-BCL2 files using the `process_srr_to_sample_symlink` function and BCL2 files using 
          the `process_bcl2_to_sample_symlink` function.
        - A shared logger is used to capture log messages from all subprocesses, and the logs are written to a
          file in the specified directory (`DATADIR`).
        - The symbolic links are created in a "results" subdirectory under the `DATADIR`, where the links follow 
          a naming pattern based on the dataset reference ID and the sample name.
    
    Example:
        main_process_srr_to_sample(args)
    """

    global DATADIR
    DATADIR = args.dir

    # === Collect all fastq.gz files
    all_files = glob.glob(f'{DATADIR}/*.fastq.gz')

    # === Split based on pattern match
    bcl2_files = [f for f in all_files if re.search(BCL2_PATTERN, f)]
    non_bcl2_files = [f for f in all_files if not re.search(BCL2_PATTERN, f)]

    # === Extract SRRs
    srr_glob = set(get_SRR_name_SRR_I_fastqgz(f) for f in non_bcl2_files)
    bcl2_glob = set(get_filename_stem(f) for f in bcl2_files)

    # === Logging setup
    shared_queue, shared_logger = get_shared_queue_and_logger_instances()
    log_file_path = get_log_file_path_symlink(DATADIR)
    shared_logger_p = multiprocessing.Process(
        target=create_shared_logger_process_for_multiprocessing,
        args=(shared_queue, log_file_path)
    )
    shared_logger_p.start()

    # debugging
    # Create a temporary logger with QueueHandler to send logs to the shared logger
    temp_logger = logging.getLogger("main_process_logger")
    temp_logger.setLevel(logging.INFO)
    temp_logger.addHandler(QueueHandler(shared_queue))
    temp_logger.propagate = False
    temp_logger.info(f"Starting temp logger")
    temp_logger.info(f"{bcl2_glob=}")

    # === Worker processes
    processes = []

    # Process non-BCL2 files
    for srr in srr_glob:
        p = multiprocessing.Process(target=process_srr_to_sample_symlink, args=(srr, DATADIR, shared_queue))
        processes.append(p)
        p.start()

    # Process BCL2-matching files separately
    for bcl2_filename in bcl2_glob:
        p = multiprocessing.Process(target=process_bcl2_to_sample_symlink, args=(bcl2_filename, DATADIR, shared_queue))
        processes.append(p)
        p.start()

    # === Wait for workers
    for p in processes:
        p.join()

    # === Shutdown logger
    shared_queue.put_nowait(None)
    shared_logger_p.join()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Setup raw data for integration model setup')
    parser.add_argument('dir', metavar='[ DATA DIR ]')
    args = parser.parse_args()

    global DATADIR
    DATADIR = args.dir

    main_process_srr_to_sample(args)
