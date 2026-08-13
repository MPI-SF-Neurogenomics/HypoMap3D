import os
import glob
import logging
import multiprocessing
import pandas as pd

from hm3d_utils.formatting import (
    File_Analysis,
    map_biosample_accession,
    extract_scx_xx,
    get_log_file_path_symlink,
    map_lane_information,
    map_sample_name,
    harmonize_bcl2_pattern
)

from hm3d_utils.logging import setup_subprocess_logger

def process_srr_to_sample_symlink(SRR: str, datadir: str, shared_queue: multiprocessing.Queue) -> None:
    """
    Processes a sample by analyzing its associated FASTQ files and creates symbolic links 
    for these files in the results directory. The symlink names follow the bcl2 standard, 
    incorporating the SCX-XX and BioSample IDs, and the lane information retrieved from 
    a metadata file. This function assumes that FASTQ files are stored in the specified 
    data directory (`datadir`) and organizes them based on sequencing read order.

    Args:
        SRR (str): The SRR (Sample Read Run) ID corresponding to the sample being processed. 
                   This is used to locate the FASTQ files associated with the sample.
        datadir (str): The directory path where the FASTQ files are stored. The function will look 
                       for files matching the SRR ID followed by `.fastq.gz` in this directory.
        shared_queue (multiprocessing.Queue): A queue from the `multiprocessing` module used for 
                                              logging or message passing, typically in parallel 
                                              processing scenarios. This queue is used to capture 
                                              logs from the subprocesses.

    Raises:
        ValueError: If the necessary FASTQ files are not found or cannot be processed.
        FileNotFoundError: If the required directories or files do not exist.
        OSError: If an error occurs during the creation of symlinks.

    Notes:
        - The function looks for FASTQ files with filenames starting with the given SRR ID 
          followed by `.fastq.gz`.
        - It sorts the files based on line length and a scoring system (t_score).
        - The symlinks are created in the 1b_bcl_FASTQ directory, which is derived from the `datadir` path 
          by replacing the last occurrence of "1a_source" with "1b_bcl_FASTQ".
        - For each sequence file, a symlink is created with a name based on the dataset reference ID, 
          BioSample ID, lane information, and read pairing (R1/R2).
        - The lane information is obtained from a metadata CSV file, where it is associated with the SRR ID.
        - If symlinks already exist, they are not overwritten.
        - The function also logs the symlink creation process using the provided `shared_queue`.

    Example:
        process_srr_to_sample_symlink("SRR12345678", "/path/to/data/directory", shared_queue)

    Notes:
        - The `datadir` argument should be the directory containing the FASTQ files associated with the given SRR ID.
        - The function automatically creates the 1b_bcl_FASTQ directory by replacing the last occurrence of "1a_source" 
          with "1b_bcl_FASTQ" in the `datadir` path.
        - Symlinks are not overwritten if they already exist.
        - The lane information for the SRR is determined by looking up the SRR in a metadata CSV file. 
          If no valid lane is found, the function defaults to lane "001".
    """
    logger = setup_subprocess_logger(shared_queue)  # Set up logger for this process

    files = glob.glob(os.path.join(datadir, f'{SRR}*.fastq.gz'))

    # Following the data directory setup, the save dir will be the 1b_bcl_FASTQ dir
    # Replace the last occurrence of "1a_source" with "1b_bcl_FASTQ" in the DATADIR path
    index = datadir.rfind("1a_source")

    if index != -1:  # Check if "1a_source" exists in the string
        save_dir = datadir[:index] + "1b_bcl_FASTQ" + datadir[index + len("1a_source"):]
    else:
        save_dir = datadir  # No replacement if "1a_source" is not found

    # Check if the directory exists and create it if not
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    # Analyze files and prepare for sorting
    filelist = [File_Analysis(f) for f in files]
    
    # Sort the files based on line length and score, descending
    sorted_filelist = sorted(filelist, key=lambda f: (f.line_length, f.t_score), reverse=True)
    
    # Map BioSample accession to SRR
    sample = map_biosample_accession(SRR, datadir) 
    # Extract dataset_ref_id from the data path
    dataset_ref_id = extract_scx_xx(datadir)
    # Map lane information
    lane = map_lane_information(SRR, datadir)

    # Create symlinks for the files in the save_dir
    for i, seq_file in enumerate(sorted_filelist):
        # Construct the target filename in the save directory (symlink name)
        if i == 0:
            symlink_name = os.path.join(save_dir, f'{dataset_ref_id}_{sample}_S1_L{lane}_R2_001.fastq.gz')
        elif i == 1:
            symlink_name = os.path.join(save_dir, f'{dataset_ref_id}_{sample}_S1_L{lane}_R1_001.fastq.gz')
        else:
            logger.info(f"{seq_file.filename} is redundant. Not creating symlink.")
            continue

        # Log the mapping from original file name to symlink
        logger.info(f"Mapping original file {seq_file.filename} to symlink {symlink_name}")
        
        # Create a symbolic link in the save directory
        if not os.path.exists(symlink_name):  # Avoid overwriting existing symlinks
            os.symlink(seq_file.filename, symlink_name)
            logger.info(f"Created symlink: {symlink_name} -> {seq_file.filename}")
        else:
            logger.warning(f"Symlink already exists: {symlink_name} -> {seq_file.filename}")



def process_bcl2_to_sample_symlink(bcl2_filename: str, datadir: str, shared_queue: multiprocessing.Queue) -> None:
    """
    Creates a symbolic link for a given BCL2 file, mapping it to a standardized sample-based
    filename in a corresponding `1b_bcl_FASTQ` directory. This function is intended to be run as a
    subprocess and logs its activity using a shared multiprocessing logging queue.

    Behavior:
        - Initializes a subprocess-safe logger using the provided shared queue.
        - Determines the full path of the input BCL2 file within the provided `datadir`.
        - Constructs the target directory path by replacing the last occurrence of "1a_source"
          in `datadir` with "1b_bcl_FASTQ` (to mirror a results directory structure).
        - Ensures the `1b_bcl_FASTQ` directory exists (creates it if needed).
        - Extracts a dataset reference ID from the `datadir` using `extract_scx_xx()`.
        - Standardizes the BCL2 filename using `harmonize_bcl2_pattern()`.
        - Maps the standardized filename to a sample name using `map_sample_name()`.
        - Constructs the full symlink path with the format: `{dataset_ref_id}_{mapped_filename}`.
        - Logs and creates the symlink in the results directory if it doesn’t already exist.

    Args:
        bcl2_filename (str): The name of the BCL2 file (e.g., "Sample01_L001_R1_001.fastq.gz").
        datadir (str): Path to the directory containing the input BCL2 file.
        shared_queue (multiprocessing.Queue): A queue for subprocess logging to avoid conflicts.

    Returns:
        None

    Notes:
        - This function avoids overwriting existing symbolic links.
        - The symlink is created in a `1b_bcl_FASTQ` directory that mirrors the input `datadir` path.
        - It is designed to be used in parallel file processing pipelines (e.g., with multiprocessing).
        - Logging is handled via a queue and consumed by a parent logging handler.

    """

    logger = setup_subprocess_logger(shared_queue)  # Set up logger for this process

    full_path_current = os.path.join(datadir, f'{bcl2_filename}')

    # Following the data directory setup, the save dir will be the 1b_bcl_FASTQ dir
    # Replace the last occurrence of "1a_source" with "1b_bcl_FASTQ" in the DATADIR path
    index = datadir.rfind("1a_source")

    if index != -1:  # Check if "1a_source" exists in the string
        save_dir = datadir[:index] + "1b_bcl_FASTQ" + datadir[index + len("1a_source"):]
    else:
        save_dir = datadir  # No replacement if "1a_source" is not found

    # Check if the directory exists and create it if not
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    #==================
    # Start creating symlink

    # Extract dataset_ref_id from the data path
    dataset_ref_id = extract_scx_xx(datadir)

    # Harmonize bcl2 pattern
    bcl2_filename = harmonize_bcl2_pattern(bcl2_filename)

    # Map BioSample/Sample name
    bcl2_filename = map_sample_name(bcl2_filename, datadir) 

    # Create symlink for the file in the save_dir
    # Construct the target filename in the save directory (symlink name)
    symlink_name = os.path.join(save_dir, f'{dataset_ref_id}_{bcl2_filename}')

    # Log the mapping from original file name to symlink
    logger.info(f"Mapping original file {full_path_current} to symlink {symlink_name}")
    
    # # Create a symbolic link in the save directory
    if not os.path.exists(symlink_name):  # Avoid overwriting existing symlinks
        os.symlink(full_path_current, symlink_name)
        logger.info(f"Created symlink: {symlink_name} -> {full_path_current}")
    else:
        logger.warning(f"Symlink already exists: {symlink_name} -> {full_path_current}")
