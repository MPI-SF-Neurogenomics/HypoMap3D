import os
import argparse
import pandas as pd

OUTPUT_FILENAME = "combined_alignment_qc_metrics_summary.csv"


# ------------------------------------------------------------------
# File discovery
# ------------------------------------------------------------------

def find_csv_files(root_dir):
    """
    Recursively search for Cell Ranger and STARsolo QC metric files.

    Cell Ranger: metrics_summary.csv
    STARsolo:    Summary.csv
    """
    target_filenames = {"metrics_summary.csv", "Summary.csv"}
    csv_files = []

    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename in target_filenames:
                csv_files.append(os.path.join(dirpath, filename))

    return csv_files


# ------------------------------------------------------------------
# Metadata inference
# ------------------------------------------------------------------

def extract_sample_name(file_path):
    """
    Extract sample name from either Cell Ranger or STARsolo paths.
    """
    parts = os.path.normpath(file_path).split(os.sep)

    if "outs" in parts:
        return parts[parts.index("outs") - 1]

    if "Solo.out" in parts:
        return parts[parts.index("Solo.out") - 1]

    return "Unknown"


def infer_alignment_technology(file_path):
    if file_path.endswith("metrics_summary.csv"):
        return "cellranger"
    if file_path.endswith("Summary.csv"):
        return "starsolo"
    return "unknown"


# ------------------------------------------------------------------
# File parsing
# ------------------------------------------------------------------

def read_cellranger_metrics(file_path):
    """
    Cell Ranger metrics_summary.csv:
    - One header row
    - One data row
    """
    return pd.read_csv(file_path)


def read_starsolo_summary(file_path):
    """
    STARsolo Summary.csv:
    - Two columns: metric | value
    - Needs to be transposed into a single-row dataframe
    """
    df = pd.read_csv(file_path, header=None)

    if df.shape[1] < 2:
        raise ValueError("STARsolo Summary.csv has fewer than two columns")

    df = df.iloc[:, :2]
    df.columns = ["metric", "value"]

    df_wide = df.set_index("metric").T
    df_wide.reset_index(drop=True, inplace=True)

    return df_wide


# ------------------------------------------------------------------
# Concatenation
# ------------------------------------------------------------------

def concatenate_csv_files(file_paths):
    dataframes = []

    for file_path in file_paths:
        try:
            alignment = infer_alignment_technology(file_path)
            sample = extract_sample_name(file_path)

            if alignment == "cellranger":
                df = read_cellranger_metrics(file_path)
            elif alignment == "starsolo":
                df = read_starsolo_summary(file_path)
            else:
                print(f"Skipping unknown file format: {file_path}")
                continue

            df["Sample"] = sample
            df["Alignment"] = alignment

            dataframes.append(df)

        except Exception as e:
            print(f"Error processing {file_path}: {e}")

    if not dataframes:
        return pd.DataFrame()

    return pd.concat(dataframes, ignore_index=True, sort=False)


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

def main(input_dir):
    csv_files = find_csv_files(input_dir)

    if not csv_files:
        print("No QC metric files found.")
        return

    combined_df = concatenate_csv_files(csv_files)

    if combined_df.empty:
        print("No valid QC data to write.")
        return

    output_path = os.path.join(input_dir, OUTPUT_FILENAME)
    combined_df.to_csv(output_path, index=False)
    print(f"Combined CSV saved to: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Combine Cell Ranger and STARsolo QC metrics into a single CSV."
    )
    parser.add_argument(
        "input_dir",
        help="Root directory to search for QC metric files."
    )

    args = parser.parse_args()
    main(args.input_dir)
