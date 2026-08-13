"""
Parallelized BAM to FASTQ conversion using Cell Ranger within a Singularity container.

This script scans a directory for BAM files and runs `cellranger bamtofastq` on each one
in parallel using Python's multiprocessing. Each BAM file is processed using a provided
Singularity image that contains the Cell Ranger binary.

Usage:
    python run_bamtofastq.py /path/to/cellranger.sif /path/to/bam_directory -j 4

Arguments:
    sif_image: Path to the Singularity image containing Cell Ranger.
    bam_dir:   Path to the directory containing BAM files to convert.
    -j / --jobs: (Optional) Number of parallel jobs to run (defaults to number of CPU cores).
"""

import argparse
import subprocess
from multiprocessing import Pool, cpu_count
from pathlib import Path


def run_bamtofastq(args):
    """
    Process a single BAM file using cellranger bamtofastq inside Singularity.

    Args:
        args (tuple): A tuple containing:
            - sif_image (Path): Path to the Singularity image.
            - bam_file (Path): Path to the BAM file.
            - output_base (Path): Directory where the FASTQ output folder will be created.
    """
    sif_image, bam_file, output_base = args
    bam_file = Path(bam_file)
    output_dir = output_base / f"{bam_file.stem}_fastq"

    if output_dir.exists():
        print(f"[SKIP] Output already exists for {bam_file.name}, skipping.")
        return

    print(f"[RUNNING] Processing {bam_file.name}...")

    try:
        subprocess.run([
            "singularity", "exec", str(sif_image),
            "cellranger", "bamtofastq",
            str(bam_file),
            str(output_dir)
        ], check=True)
        print(f"[DONE] Finished {bam_file.name}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed processing {bam_file.name}: {e}")


def parse_args():
    """
    Parse command-line arguments using argparse.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Parallel Cell Ranger bamtofastq using Singularity"
    )
    parser.add_argument(
        "sif_image", type=Path,
        help="Path to the Singularity image containing Cell Ranger"
    )
    parser.add_argument(
        "bam_dir", type=Path,
        help="Directory containing BAM files"
    )
    parser.add_argument(
        "-j", "--jobs", type=int, default=cpu_count(),
        help="Number of parallel jobs (default: number of CPU cores)"
    )
    return parser.parse_args()


def main():
    """
    Main function that manages the parallel execution of bamtofastq
    across all BAM files in the specified directory.
    """
    args = parse_args()

    if not args.sif_image.is_file():
        print(f"Error: Singularity image not found at {args.sif_image}")
        return

    if not args.bam_dir.is_dir():
        print(f"Error: BAM directory not found at {args.bam_dir}")
        return

    bam_files = list(args.bam_dir.glob("*.bam"))

    if not bam_files:
        print(f"No BAM files found in {args.bam_dir}")
        return

    print(f"Found {len(bam_files)} BAM files. Running with {args.jobs} workers...")

    arg_list = [(args.sif_image, bam_file, args.bam_dir) for bam_file in bam_files]

    with Pool(processes=min(args.jobs, len(bam_files))) as pool:
        pool.map(run_bamtofastq, arg_list)

    print("All jobs completed.")


if __name__ == "__main__":
    main()
