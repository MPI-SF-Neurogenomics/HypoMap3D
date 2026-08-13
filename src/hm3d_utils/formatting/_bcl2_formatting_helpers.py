
import os 
import re 
import glob
import pandas as pd
from pandas.errors import ParserError
import csv

def get_log_file_path_symlink(datadir: str) -> str:
    """
    Derives the log file path by replacing '1a_source' with '1b_bcl_FASTQ/logs' in the DATADIR path.
    Prepends the SCX-XX string to the log file name.
    
    Returns:
        str: The path where the log file will be saved.
    """
    
    # Extract the SCX-XX information from the DATADIR path
    dataset_ref_id = extract_scx_xx(datadir)
    
    # Find the index of "data" in DATADIR
    index = datadir.rfind("1a_source")

    if index != -1:  # Check if "1a_source" exists in the string
        log_dir = datadir[:index] + "1b_bcl_FASTQ/logs"
    else:
        log_dir = datadir  # No replacement if "1a_source" is not found
    
    # Ensure the log directory exists
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    
    # Derive the log file path with SCX-XX prepended
    log_file_path = os.path.join(log_dir, f"{dataset_ref_id}_symlink_creation.log")
    return log_file_path


def get_SRR_name_SRR_I_fastqgz(file_path: str) -> str:
    """
    Extracts the SRR (Sequence Read Run) ID from a given FASTQ file path.

    This function assumes that the filename contains the SRR ID as the prefix, followed by an underscore and read information.

    Args:
        file_path (str): Path to the FASTQ file.

    Returns:
        str: The extracted SRR ID, e.g., "SRR24442422".

    Raises:
        ValueError: If the SRR ID cannot be found in the filename.
    """
    base = os.path.basename(file_path)
    filename_without_extension = base.split('.')[0]
    
    # Use regular expression to capture the SRR number (with SRR prefix and number)
    match = re.match(r'(SRR\d+)', filename_without_extension)
    
    if match:
        return match.group(1)  # Return the SRR ID (e.g., SRR24442422)
    else:
        raise ValueError(f"No Matching SRR ID found in {filename_without_extension}")

    
def extract_prefix_from_filename(file_path: str, pattern: str) -> str:
    """
    Extracts a user-defined pattern from a given filename.

    Args:
        file_path (str): Path to the FASTQ file.
        pattern (str): Regex pattern to extract the desired part of the filename.

    Returns:
        str: The extracted prefix.

    Raises:
        ValueError: If the pattern does not match.
    """
    base = os.path.basename(file_path)
    filename_without_extension = base.split('.')[0]
    
    match = re.match(pattern, filename_without_extension)
    if match:
        return match.group(1)
    else:
        raise ValueError(f"No match found for pattern '{pattern}' in '{filename_without_extension}'")

def map_biosample_accession(SRR: str, datadir: str) -> str:
    """
    Maps the given SRR (Sequence Read Run) ID to its corresponding BioSample accession by searching
    for a metadata file in the specified directory (datadir) and reading it.

    Args:
        SRR (str): The SRR ID to map. This is the Sequence Read Run ID, typically in the format 'SRRXXXXX'.
        datadir (str): The directory path where the metadata file is located. 
                       This directory is searched for a CSV file containing metadata.

    Returns:
        str: The BioSample accession associated with the given SRR. This value is extracted from 
             the metadata file.

    Raises:
        FileNotFoundError: If no metadata file matching "SraRunTableMetadata" (case-insensitive) is 
                           found in the specified `datadir`.
        KeyError: If the SRR is not found in the metadata file. This could happen if the SRR ID 
                  is not listed in the 'Run' column of the metadata file.

    Example:
        biosample = map_biosample_accession("SRR12345678", "/path/to/metadata/directory")
    
    Notes:
        - This function performs a case-insensitive search for files containing "SraRunTableMetadata" 
          in their name (e.g., "SraRunTableMetadata.csv", "sraruntablemetadata_v2.CSV", etc.).
        - The metadata file should contain at least two columns: 'Run' (for SRR IDs) and 
          'BioSample' (for BioSample accessions).
        - If multiple metadata files match, the first one encountered is used.
    """
    # Use glob to list all files in the directory
    all_files = glob.glob(os.path.join(datadir, "*"))

    # Regex pattern to match files containing 'SraRunTableMetadata' regardless of case
    pattern = re.compile(r"sraruntablemetadata", re.IGNORECASE)

    # Filter files matching the pattern
    metadata_files = [f for f in all_files if pattern.search(os.path.basename(f))]

    if not metadata_files:
        raise FileNotFoundError(f"No metadata file found in {datadir} that matches 'SraRunTableMetadata'.")

    # Use the first match
    metadata_file = metadata_files[0]

    try:
        # read csv
        try:
            biosample_mapper = pd.read_csv(metadata_file, index_col="Run")["BioSample"]
        except (ParserError, ValueError):
            biosample_mapper = pd.read_csv(metadata_file, index_col="Run", sep=';')["BioSample"]
        return biosample_mapper[SRR]
    except KeyError:
        raise KeyError(f"SRR {SRR} not found in the metadata file.")

def map_lane_information(SRR: str, datadir: str) -> str:
    """
    Maps a Sequence Read Run (SRR) ID to a lane number using sample metadata from a supported file
    within the specified directory.

    Supported metadata formats:
        - **SRA Run Table Metadata** (CSV files containing 'SraRunTableMetadata' in the filename)
        - **NEMO Manifest** (CSV/TSV files containing 'nemo_manifest' in the filename)

    Behavior:
        1. The function loads the metadata file using `load_metadata_file()`, which returns both
           the DataFrame and a string indicating the source format.
        2. Depending on the source:
            - If `"sraruntablemetadata"`: the DataFrame is indexed by the 'Run' column.
            - If `"nemo"`: the DataFrame is indexed by the 'sample_id' column.
        3. If a 'lane' or 'Lane' column is present, its value for the given SRR is used as the lane number.
            - It is coerced to an integer, checked for validity (1–999), and returned as a
              zero-padded string (e.g., "001").
            - If the value is invalid, missing, or out of bounds, the function returns "001".
        4. If no explicit lane column exists:
            - Lane is computed dynamically:
                - The BioSample associated with the SRR is identified.
                - All SRRs with the same BioSample are collected and sorted by the numeric
                  suffix of their SRR IDs.
                - The SRR's index in this sorted list determines its lane number (starting from "001").
        5. If the required 'BioSample' column is missing during dynamic computation, a KeyError is raised.

    Args:
        SRR (str): The Sequence Read Run ID to map (e.g., "SRR26102878").
        datadir (str): Path to the directory containing the metadata file.

    Returns:
        str: A three-digit zero-padded string representing the lane number for the given SRR.

    Raises:
        FileNotFoundError: If no valid metadata file is found in the specified directory.
        KeyError:
            - If expected columns ('Run', 'sample_id', or 'BioSample') are missing.
            - If the specified SRR is not present in the metadata file.
        RuntimeError: If the metadata format is not recognized.

    Example:
        >>> map_lane_information("SRR26102878", "/path/to/metadata")
        '001'

    Notes:
        - SRRs are sorted numerically based on the digits in their ID (e.g., "SRR10" > "SRR2").
        - If an SRR is the only one for its BioSample, it is always assigned lane "001".
        - In the absence of explicit lane information, dynamic assignment ensures SRRs from the
          same BioSample get consecutive lane numbers starting from "001".

    Related Functions:
        - load_metadata_file(datadir): Loads metadata and detects its source type.
        - try_read_file(filepath): Robustly reads CSV/TSV files with auto delimiter detection.
    """


    df, source = load_metadata_file(datadir)

    if source == "sraruntablemetadata":
        # Make sure 'Run' is a column
        if 'Run' not in df.columns:
            raise KeyError("'Run' column is missing from the metadata file.")

        df.set_index('Run', inplace=True)

    elif source == "nemo":
        if 'sample_id' not in df.columns:
            raise KeyError("'sample_id' column is missing from the metadata file.")

        df.set_index('sample_id', inplace=True)

    else:
        raise RuntimeError("Unexpected metadata type.")

    # Check for explicit lane column
    lane_column = None
    if "lane" in df.columns:
        lane_column = "lane"
    elif "Lane" in df.columns:
        lane_column = "Lane"

    if lane_column:
        try:
            lane = int(df.loc[SRR, lane_column])
            if not (1 <= lane <= 999):
                lane = 1
            return f"{lane:03}"
        except (KeyError, ValueError, TypeError):
            return "001"

    # No lane column, compute dynamically
    if 'BioSample' not in df.columns:
        raise KeyError("'BioSample' column is missing for dynamic lane assignment.")

    # Get BioSample for given SRR
    try:
        biosample = df.loc[SRR, 'BioSample']
    except KeyError:
        raise KeyError(f"SRR {SRR} not found in metadata.")

    # Get all Run IDs for that BioSample
    same_sample_df = df[df['BioSample'] == biosample].copy()

    # Extract numeric part of SRRs to sort properly
    same_sample_df['Run_numeric'] = same_sample_df.index.str.extract(r'(\d+)$').astype(int)

    # Sort by numeric Run ID
    same_sample_df.sort_values('Run_numeric', inplace=True)

    # Re-index from 1 and assign lanes
    run_to_lane = {run: idx + 1 for idx, run in enumerate(same_sample_df.index)}

    # Lookup lane for the target SRR
    lane = run_to_lane.get(SRR, 1)
    return f"{lane:03}"




def extract_scx_xx(file_path: str) -> str:
    """
    Extracts the SCX-XX string (e.g., "SC2-23") from a given file path. The SCX-XX format is assumed to be
    part of the directory name, and the regular expression looks for this pattern.

    Args:
        file_path (str): The file path containing the SCX-XX string.

    Returns:
        str: The extracted SCX-XX string.

    Raises:
        ValueError: If the SCX-XX string is not found in the file path.
    """
    # Regular expression to match 'SC' followed by a number, a hyphen, and another number
    match = re.search(r'SC\d+-\d+', file_path)
    
    if match:
        return match.group(0)  # Return the matched SCX-XX part
    else:
        # Raise an error if SCX-XX information is not found in the file path
        raise ValueError("No SCX-XX information in DATADIR variable passed. Check the structure of your directory: The data should be in a directory named after its ref ID")


def harmonize_bcl2_pattern(file_path: str) -> str:
    """
    Harmonizes a BCL2 FASTQ filename pattern by ensuring all components are separated with underscores.

    Converts things like:
        "NW_TX0057-6_S01_L003.I1_001.fastq.gz" -> "NW_TX0057-6_S01_L003_I1_001.fastq.gz"
        "NW_TX0057-6_L003.I1_001.fastq.gz"    -> "NW_TX0057-6_S01_L003_I1_001.fastq.gz"

    Args:
        file_path (str): The path or filename to harmonize.

    Returns:
        str: The harmonized filename (not full path).
    """
    # Main pattern with S number
    full_pattern = re.compile(r'(S\d+)[_\-\.](L\d+)[_\-\.]([RI]\d+)[_\-\.](\d+)(\.fastq\.gz)')
    match = full_pattern.search(file_path)

    if match:
        new_bcl2_part = f"{match.group(1)}_{match.group(2)}_{match.group(3)}_{match.group(4)}{match.group(5)}"
        harmonized_filename = full_pattern.sub(new_bcl2_part, file_path)
        return harmonized_filename

    # Fallback pattern without S number
    fallback_pattern = re.compile(r'(L\d+)[_\-\.]([RI]\d+)[_\-\.](\d+)(\.fastq\.gz)')
    fallback_match = fallback_pattern.search(file_path)

    if fallback_match:
        # Inject default S01
        new_bcl2_part = f"S01_{fallback_match.group(1)}_{fallback_match.group(2)}_{fallback_match.group(3)}{fallback_match.group(4)}"
        harmonized_filename = fallback_pattern.sub(new_bcl2_part, file_path)
        return harmonized_filename

    # No matches found
    return file_path


def get_filename_stem(filepath: str) -> str:
    """
    Extracts the filename without directory path.
    
    Example:
        Input:  "/data/sample1_S1_L001_R1_001.fastq.gz"
        Output: "sample1_S1_L001_R1_001.fastq.gz"
    """
    return os.path.basename(filepath)


def map_sample_name(filename: str, datadir: str) -> str:
    """
    Maps the sample name prefix of a harmonized BCL2 filename to a new name using a CSV-based mapping.

    This function identifies the BCL2 pattern in a given harmonized filename, extracts the sample prefix 
    (i.e., everything before the BCL2 pattern), maps it to a new sample name using the metadata CSV file 
    from the specified directory, and then reconstructs the filename with the new sample name prefix 
    while keeping the original BCL2 pattern suffix intact.

    Args:
        filename (str): Harmonized BCL2 FASTQ filename (not the full path). The filename must match the BCL2 pattern.
        datadir (str): Directory path where the CSV mapping file is located. This directory should contain metadata files 
                       to map the sample prefix to a new name.

    Returns:
        str: New filename with the mapped sample name prefix and the original BCL2 pattern suffix.

    Raises:
        ValueError: If the filename does not match the expected BCL2 pattern, or if the mapping fails.
        FileNotFoundError: If no valid metadata file is found in the directory.
        RuntimeError: If an unexpected metadata type is encountered during the mapping process.
    """
    # BCL2 pattern to identify where the pattern part starts
    bcl2_pattern = re.compile(r'(_S\d+_L\d+_[RI]\d+_\d+\.fastq\.gz)$')

    # Try to find the BCL2 suffix (e.g., _S01_L003_I1_001.fastq.gz)
    match = bcl2_pattern.search(filename)
    if not match:
        raise ValueError(f"Filename '{filename}' does not match expected BCL2 pattern.")

    # Extract the BCL2 pattern suffix and the original prefix (everything before it)
    bcl2_suffix = match.group(1)
    original_prefix = filename[:match.start(1)]

    # Map the original prefix to a new sample name using the appropriate metadata
    mapped_prefix = get_sample_name(datadir, original_prefix)

    # Reconstruct the final filename by combining the mapped prefix with the original BCL2 suffix
    mapped_filename = f"{mapped_prefix}{bcl2_suffix}"
    
    return mapped_filename


def try_read_file(filepath: str) -> pd.DataFrame:
    """
    Attempts to read a CSV file, automatically detecting the delimiter 
    or using a fallback delimiter if detection fails.

    This function tries to read the given file using `pandas.read_csv()`. 
    It first attempts to detect the separator using `csv.Sniffer`. 
    If detection fails, it retries with a semicolon-separated format. 
    If both attempts fail, it raises a `ValueError` indicating the failure.

    Args:
        filepath (str): The path to the CSV file to read.

    Returns:
        pd.DataFrame: The DataFrame read from the file.

    Raises:
        ValueError: If both attempts to read the file fail (auto-detection, semicolon).
    """
    try:
        # Try to automatically detect the separator using csv.Sniffer
        with open(filepath, 'r', encoding='utf-8') as file:
            sample = file.read(1024)  # Read a sample of the file (first 1024 bytes)
            delimiter = csv.Sniffer().sniff(sample).delimiter  # Detect delimiter

        # Try reading the CSV with the detected separator
        return pd.read_csv(filepath, sep=delimiter)
    
    except Exception as e:
        # If detection fails or any other error, try semicolon separator
        try:
            return pd.read_csv(filepath, sep=';')
        except Exception as e:
            # If semicolon also fails, try default (comma) separator
            try:
                return pd.read_csv(filepath)
            except Exception as e:
                raise ValueError(f"Failed to read metadata file: {filepath}") from e

def load_metadata_file(datadir: str) -> tuple[pd.DataFrame, str]:
    """
    Loads a metadata file from the given directory and determines its source type.

    Returns:
        Tuple of (DataFrame, source_type), where source_type is either 'sraruntablemetadata' or 'nemo'.

    Raises:
        FileNotFoundError: If no appropriate metadata file is found.
    """
    for filename in os.listdir(datadir):
        if "sraruntablemetadata" in filename.lower() and filename.lower().endswith(".csv"):
            return try_read_file(os.path.join(datadir, filename)), "sraruntablemetadata"

    for filename in os.listdir(datadir):
        name_lower = filename.lower()
        if "nemo_manifest" in name_lower and (name_lower.endswith(".csv") or name_lower.endswith(".tsv")):
            path = os.path.join(datadir, filename)
            if name_lower.endswith(".tsv"):
                return pd.read_csv(path, sep="\t"), "nemo"
            else:
                return try_read_file(path), "nemo"

    raise FileNotFoundError("No valid metadata file found in directory.")




def find_column_case_insensitive(columns, target):
    """
    Searches for a column name in a list of columns, ignoring case.

    This function checks if the target string (column name) exists in the 
    given columns list, case-insensitively. It will return the first match it finds 
    and ignore case differences.

    Args:
        columns (list): A list of column names (strings) to search through.
        target (str): The target column name to find, case-insensitively.

    Returns:
        str: The column name that matches the target (case-insensitive match).

    Raises:
        ValueError: If no match is found.
    """
    for col in columns:
        if col.lower() == target.lower():
            return col
    return None



def map_prefix_from_nemomanifest(df: pd.DataFrame, original_prefix: str) -> str:
    """
    Maps a prefix to its sample ID using a NeMO manifest DataFrame.

    Args:
        df (pd.DataFrame): The loaded NeMO manifest.
        original_prefix (str): The prefix string to search for (as substring).

    Returns:
        str: The corresponding sample_id value.

    Raises:
        ValueError: If mapping cannot be performed.
    """
    # Potential columns to use as index, in priority order
    candidate_index_names = ["file_id", "file_name", "filename"]

    index_col = None
    for name in candidate_index_names:
        index_col = find_column_case_insensitive(df.columns, name)
        if index_col:
            break

    if not index_col:
        raise ValueError("Could not find a valid index column ('file_id', 'file_name', or 'filename') in NeMO manifest.")

    df.set_index(index_col, inplace=True)

    # Find the first row in the index where the original_prefix is a substring
    matching_indices = [
        idx for idx in df.index
        if original_prefix in str(idx)
    ]

    if not matching_indices:
        raise ValueError(f"No entry in column '{index_col}' contains original prefix '{original_prefix}'")

    # Now find the sample_id column
    sample_id_col = find_column_case_insensitive(df.columns, "sample_id")
    if not sample_id_col:
        raise ValueError("Could not find 'sample_id' column in NeMO manifest.")

    mapped_value = df.loc[matching_indices[0], sample_id_col]
    return str(mapped_value)



def get_sample_name(datadir: str, original_prefix: str) -> str:
    """
    Maps an original prefix to a sample name using metadata from a CSV file.

    This function loads the appropriate metadata file from the provided directory
    (`datadir`), determines the correct metadata source (either `sraruntablemetadata` or `nemo`), 
    and maps the `original_prefix` to a corresponding sample name (`BioSample` or `sample_id`) 
    from the metadata file.

    Args:
        datadir (str): Directory containing the metadata file.
        original_prefix (str): The prefix extracted from the filename to map to a sample.

    Returns:
        str: The mapped sample name (e.g., BioSample or sample_id).

    Raises:
        FileNotFoundError: If no valid metadata file is found in the directory.
        ValueError: If the original prefix cannot be mapped to a sample name.
        RuntimeError: If an unexpected metadata type is encountered.
    """
    df, source = load_metadata_file(datadir)

    if source == "sraruntablemetadata":
        return map_prefix_from_sraruntablemetadata(df, original_prefix)

    elif source == "nemo":
        return map_prefix_from_nemomanifest(df, original_prefix)

    raise RuntimeError("Unexpected metadata type.")

def map_prefix_from_sraruntablemetadata(df, original_prefix: str) -> str:
    """
    Maps an original prefix to a new sample name using the SRARunTableMetadata CSV file.

    This function first searches for a column containing the `original_prefix` in the dataframe, 
    then sets the index to that column. It then looks for the "BioSample" column (case-insensitive)
    and uses the `original_prefix` to retrieve the corresponding sample name from that column.

    Args:
        df (pd.DataFrame): The dataframe containing SRARunTableMetadata.
        original_prefix (str): The prefix to map, typically extracted from the filename.

    Returns:
        str: The mapped sample name corresponding to the original prefix.

    Raises:
        ValueError: If the original prefix cannot be found in the relevant column.
        ValueError: If the "BioSample" column cannot be found.
        ValueError: If the original prefix is not found in the index column.
    """
    # Find the column containing original_prefix
    matching_columns = []
    for col in df.columns:
        if df[col].astype(str).str.strip().eq(original_prefix).any():
            matching_columns.append(col)

    if not matching_columns:
        raise ValueError(f"Could not find a column containing the original prefix: '{original_prefix}'")

    index_col = matching_columns[0]
    df.set_index(index_col, inplace=True)

    # Find the "BioSample" column (case-insensitive)
    biosample_col = find_column_case_insensitive(df.columns, "BioSample")
    if not biosample_col:
        raise ValueError("BioSample column not found in SRARunTableMetadata file.")

    try:
        # Map the original prefix to the BioSample
        mapped_value = df.loc[original_prefix, biosample_col]
        return str(mapped_value)
    except KeyError:
        raise ValueError(f"Original prefix '{original_prefix}' not found in index column '{index_col}'.")

