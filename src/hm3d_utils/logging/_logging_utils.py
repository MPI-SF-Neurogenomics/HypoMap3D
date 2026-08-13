
import os
import re

def get_log_file_path_symlink_forwarding(results_dir: str) -> str:
    """
    Derives the log file path by replacing 'results' with 'logs' in the results_dir path.
    
    Returns:
        str: The path where the log file will be saved.
    """
    
    # Find the index of "data" in DATADIR
    index = results_dir.rfind("results")

    if index != -1:  # Check if "data" exists in the string
        log_dir = results_dir[:index] + "logs" + results_dir[index + len("results"):]
    else:
        log_dir = results_dir  # No replacement if "data" is not found
    
    # Ensure the log directory exists
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    
    # Derive the log file path with SCX-XX prepended
    log_file_path = os.path.join(log_dir, f"symlink_forwarding.log")
    return log_file_path