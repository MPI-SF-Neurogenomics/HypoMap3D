#!/usr/bin/env python3

import os
import shutil
import argparse

def get_unique_filename(target_dir, filename):
    """
    Generate a unique filename in the target directory by appending a numeric suffix
    if the file already exists.

    Args:
        target_dir (str): Directory to check for conflicts.
        filename (str): Desired filename.

    Returns:
        str: A filename guaranteed not to conflict.
    """
    base, ext = os.path.splitext(filename)
    counter = 1
    new_filename = filename

    while os.path.exists(os.path.join(target_dir, new_filename)):
        new_filename = f"{base}_{counter}{ext}"
        counter += 1

    return new_filename

def move_files_to_target(target_dir: str, prepend_subdir: bool = False, dry_run: bool = False, delete_empty_subdirs: bool = False):
    """
    Moves all files from the subdirectories of the target directory into the target directory itself.
    
    Optionally, removes empty subdirectories after moving files.

    Args:
        target_dir (str): The directory where all files will be collected.
        prepend_subdir (bool): Whether to prepend subdirectory names to filenames.
        dry_run (bool): If True, simulates actions without modifying the filesystem.
        delete_empty_subdirs (bool): If True, deletes empty subdirectories after moving files.
    """
    target_dir = os.path.abspath(target_dir)
    print(f"Target directory: {target_dir}")
    if dry_run:
        print("Dry run mode: No files will be moved.\n")

    # Simulated filename tracker
    future_filenames = set(os.listdir(target_dir)) if dry_run else None

    for root, dirs, files in os.walk(target_dir):
        if root == target_dir:
            continue

        subdir_name = os.path.basename(root)

        for original_filename in files:
            old_path = os.path.join(root, original_filename)

            # Work on a copy of the filename (do not modify original_filename)
            filename = f"{subdir_name}_{original_filename}" if prepend_subdir else original_filename

            # Determine unique filename
            base, ext = os.path.splitext(filename)
            counter = 1
            new_filename = filename

            while True:
                if dry_run:
                    if new_filename not in future_filenames:
                        break
                else:
                    if not os.path.exists(os.path.join(target_dir, new_filename)):
                        break

                new_filename = f"{base}_{counter}{ext}"
                counter += 1

            new_path = os.path.join(target_dir, new_filename)

            if dry_run:
                print(f"[Dry run] Would move: {old_path} → {new_path}")
                future_filenames.add(new_filename)  # Simulate the move
            else:
                print(f"Moving: {old_path} → {new_path}")
                shutil.move(old_path, new_path)

        # Optionally remove empty subdirectories after files have been moved
        if delete_empty_subdirs and not os.listdir(root):  # If directory is empty
            print(f"Removing empty subdirectory: {root}")
            os.rmdir(root)

    if dry_run:
        print("\nDry run complete. No files were moved.")
    else:
        print("\nDone moving files.")

if __name__ == "__main__":
    """
    Script to flatten a directory by moving all files from its subdirectories into the top-level directory.

    Usage:
        python move_subdir_files_into_target_dir.py /path/to/target [--prepend-subdir] [--dry-run] [--delete-empty-subdirs]

    Positional Arguments:
        target              Path to the target directory to consolidate files into.

    Optional Arguments:
        --prepend-subdir    Prepend subdirectory name to each filename to avoid name collisions.
        --dry-run           Simulate file operations without making changes.
        --delete-empty-subdirs  Remove empty subdirectories after moving files.

    Example:
        python move_subdir_files_into_target_dir.py ./data --prepend-subdir --dry-run --delete-empty-subdirs
    """
    parser = argparse.ArgumentParser(description="Move all files from subdirectories into a target directory.")
    parser.add_argument("target", type=str, help="Path to the target directory.")
    parser.add_argument("--prepend-subdir", action="store_true",
                        help="Prepend the subdirectory name to each filename.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Simulate file operations without making changes.")
    parser.add_argument("--delete-empty-subdirs", action="store_true",
                        help="Delete empty subdirectories after moving files.")
    
    args = parser.parse_args()
    move_files_to_target(args.target, args.prepend_subdir, args.dry_run, args.delete_empty_subdirs)
