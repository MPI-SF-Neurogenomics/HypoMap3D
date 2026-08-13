#!/usr/bin/env python3

import os
import tarfile
import argparse
from concurrent.futures import ThreadPoolExecutor

def extract_tar(tar_file, extract_dir):
    """
    Extracts a tar file to the specified directory.

    This function opens a tar file, extracts its contents, and places them in the provided
    extraction directory.

    Args:
        tar_file (str): Path to the tar file.
        extract_dir (str): Directory to extract the contents to.

    Raises:
        Exception: If an error occurs while extracting the tar file.
    """
    try:
        with tarfile.open(tar_file, 'r') as tar:
            print(f"Extracting: {tar_file}")
            tar.extractall(path=extract_dir)
    except Exception as e:
        print(f"Error extracting {tar_file}: {e}")


def find_fastq_tar_files(directory):
    """
    Finds all .fastq.tar files in the specified directory.

    This function searches through the directory (including subdirectories) to find all files
    that end with the '.fastq.tar' extension.

    Args:
        directory (str): Path to the directory to search in.

    Returns:
        list: A list of file paths for all .fastq.tar files found.
    """
    tar_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".fastq.tar"):
                tar_files.append(os.path.join(root, file))
    return tar_files


def main(directory):
    """
    Main function to find and extract .fastq.tar files in the given directory in parallel.

    This function first checks if the provided directory exists. Then, it searches for all
    .fastq.tar files within that directory (and subdirectories). It uses threading to parallelize
    the extraction of those files, thus speeding up the extraction process.

    Args:
        directory (str): Path to the directory to search in and extract files.

    Raises:
        FileNotFoundError: If the directory does not exist.
        ValueError: If no .fastq.tar files are found in the directory.
    """
    if not os.path.isdir(directory):
        print(f"Error: Directory {directory} does not exist.")
        return

    tar_files = find_fastq_tar_files(directory)

    if not tar_files:
        print(f"No .fastq.tar files found in {directory}.")
        return

    # Use ThreadPoolExecutor to parallelize extraction
    with ThreadPoolExecutor() as executor:
        # Submit the extraction tasks to the executor
        for tar_file in tar_files:
            executor.submit(extract_tar, tar_file, directory)

    print("Extraction complete!")


if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Extract .fastq.tar files in a directory in parallel.")
    parser.add_argument("directory", type=str, help="Directory containing the .fastq.tar files to extract.")
    args = parser.parse_args()

    # Run the main function
    main(args.directory)
