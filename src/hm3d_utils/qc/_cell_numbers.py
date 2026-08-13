
import pandas as pd
from typing import List

import seaborn as sns
import matplotlib.pyplot as plt

def summarize_qc_counts(
    doublet_file: str,
    qc_dir: str,
    datasets: List[str],
) -> pd.DataFrame:
    """
    Summarizes quality control and doublet filtering information for a set of datasets.

    Parameters:
        doublet_file (str): Path to the CSV file containing doublet predictions.
        qc_dir (str): Path to the directory containing QC summary CSVs.
        datasets (List[str]): List of dataset identifiers (strings).

    Returns:
        pd.DataFrame: A summary DataFrame with QC and doublet counts per sample.
    """

    # Load doublet results
    doublets = pd.read_csv(doublet_file)
    doublets["Sample"] = doublets["Unnamed: 0"].apply(lambda x: "_".join(x.split("_")[1:]))

    # Summarize doublets
    doublet_qc = doublets.groupby("Sample")["DoubletToRemove"].agg(
        num_cells_pre_doublet_removal="count",
        num_doublets_removed="sum"
    )

    # Aggregate QC metrics across datasets
    qc_summaries = []

    for ds in datasets:
        qc_path = f"{qc_dir}/{ds}_qc_summary.csv"
        qc_df = pd.read_csv(qc_path)

        failed_qc = qc_df.groupby("Sample")["failed_qc"].agg(
            num_cells_pre_qc="count",
            num_cells_removed_via_qc_threshold="sum"
        )

        manually_removed = qc_df.groupby("Sample")["manually_removed"].agg(
            num_cells_removed_manually="sum"
        )

        # Combine QC summaries for this dataset
        dataset_qc = pd.concat([failed_qc, manually_removed], axis=1)
        qc_summaries.append(dataset_qc)

    # Combine all datasets' QC summaries
    qc_all = pd.concat(qc_summaries, axis=0)

    # Merge QC and doublet summaries
    final_qc = pd.concat([doublet_qc, qc_all], axis=1)

    # Compute number of cells retained after doublet removal
    final_qc["num_cells_retained"] = (
        final_qc["num_cells_pre_doublet_removal"] - final_qc["num_doublets_removed"]
    )

    # Reorder columns
    final_qc = final_qc.loc[:, [
        "num_cells_pre_qc",
        "num_cells_removed_via_qc_threshold",
        "num_cells_removed_manually",
        "num_cells_pre_doublet_removal",
        "num_doublets_removed",
        "num_cells_retained"
    ]]

    return final_qc

def melt_qc_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Reshapes a QC summary DataFrame from wide to long format for visualization.

    This function:
      - Selects key QC stage columns (`num_cells_pre_qc`, `num_cells_pre_doublet_removal`, `num_cells_retained`)
      - Renames them to more intuitive labels (`initial`, `post_qc`, `post_doublet`)
      - Melts the DataFrame into long format with columns: "Sample", "qc_stage", "num_cells"
      - Extracts the dataset name from each sample ID (assuming samples are named like "Dataset_Sample")

    Parameters:
        df (pd.DataFrame): A QC summary DataFrame with "Sample" as index and the columns:
                           "num_cells_pre_qc", "num_cells_pre_doublet_removal", "num_cells_retained".

    Returns:
        pd.DataFrame: A melted DataFrame with columns:
                      ["Sample", "qc_stage", "num_cells", "Dataset"],
                      suitable for plotting with seaborn/matplotlib.
    """
    qc_long = df.loc[:, ["num_cells_pre_qc", "num_cells_pre_doublet_removal", "num_cells_retained"]].reset_index()

    qc_long.columns = ["Sample", "initial", "post_qc", "post_doublet"]

    qc_long_melted = pd.melt(
        qc_long,
        id_vars="Sample",
        var_name="qc_stage",
        value_name="num_cells"
    )

    qc_long_melted["Dataset"] = qc_long_melted["Sample"].apply(lambda x: "_".join(x.split("_")[:-1]))

    return qc_long_melted



def plot_per_sample_qc_barplot(
    df,
    save_path=None,
    x="Sample",
    y="num_cells",
    hue="qc_stage",
    title="Number of cells per sample throughout the pipeline",
    height=6,
    aspect=2.5,
    dpi=300
):
    """
    Create and save a rasterized Seaborn barplot showing QC stages per sample.

    Parameters:
        df (pd.DataFrame): Melted DataFrame containing long-format QC data.
        save_path (str or None): If provided, saves the figure as a rasterized PDF at this path.
        x (str): Column name to use for x-axis (default: "Sample").
        y (str): Column name to use for y-axis (default: "num_cells").
        hue (str): Column name to use for hue grouping (default: "qc_stage").
        title (str): Title of the plot.
        height (float): Height of the plot.
        aspect (float): Aspect ratio of the plot (width = height × aspect).
        dpi (int): Dots per inch for saved figure.
    """

    # Create the catplot
    g = sns.catplot(
        data=df,
        x=x,
        y=y,
        hue=hue,
        kind="bar",
        height=height,
        aspect=aspect
    )

    # Set plot title
    g.fig.suptitle(title, fontsize=16)

    # Adjust spacing for title
    g.fig.subplots_adjust(top=0.9)

    # Rotate x-axis labels
    g.set_xticklabels(rotation=90, horizontalalignment='right')

    # Move legend outside
    g._legend.set_bbox_to_anchor((1.1, 0.8))
    g._legend.set_title(hue)

    # Add annotations inside bars
    for ax in g.axes.flat:
        for p in ax.patches:
            height = p.get_height()
            if height > 0:
                ax.annotate(
                    f"{int(height)}",
                    (p.get_x() + p.get_width() / 2, height / 2),
                    ha='center',
                    va='center',
                    fontsize=8,
                    rotation=90,
                    color='white'
                )

    # Rasterize the bars for smaller PDF size
    for ax in g.axes.flat:
        for bar in ax.patches:
            bar.set_rasterized(True)

    # Save if requested
    if save_path is not None:
        plt.tight_layout()
        g.savefig(save_path, dpi=dpi, bbox_inches="tight", format="pdf")
        
    plt.show()




def plot_per_dataset_qc_barplot(
    df,
    save_path=None,
    title="Number of cells per dataset throughout the pipeline",
    hue_order=("initial", "post_qc", "post_doublet"),
    height=6,
    aspect=1.5,
    dpi=300
):
    """
    Create and optionally save a rasterized Seaborn barplot summarizing QC stages per dataset.

    Parameters:
        df (pd.DataFrame): DataFrame with columns ["Dataset", "qc_stage", "num_cells"] in long format.
        save_path (str or None): If provided, saves the figure as a rasterized PDF at this path.
        title (str): Title for the plot.
        hue_order (tuple): Order of QC stages in the legend.
        height (float): Height of the plot in inches.
        aspect (float): Aspect ratio of the plot (width = height × aspect).
        dpi (int): Resolution (DPI) for saved figure.
    """

    # Aggregate total number of cells per dataset and QC stage
    df_agg = df[["Dataset", "qc_stage", "num_cells"]].groupby(["Dataset", "qc_stage"]).sum().reset_index()

    # Create the catplot
    g = sns.catplot(
        data=df_agg,
        x="Dataset",
        y="num_cells",
        hue="qc_stage",
        hue_order=hue_order,
        kind="bar",
        height=height,
        aspect=aspect
    )

    # Set plot title
    g.fig.suptitle(title, fontsize=16)
    g.fig.subplots_adjust(top=0.9)

    # Rotate x-axis labels
    g.set_xticklabels(rotation=90, horizontalalignment='right')

    # Move legend outside
    g._legend.set_bbox_to_anchor((1.1, 0.8))
    g._legend.set_title("QC Stage")

    # Add annotations inside bars
    for ax in g.axes.flat:
        for p in ax.patches:
            height = p.get_height()
            if height > 0:
                ax.annotate(
                    f"{int(height)}",
                    (p.get_x() + p.get_width() / 2, height / 2),
                    ha='center',
                    va='center',
                    fontsize=11,
                    rotation=90,
                    color='white'
                )

    # Rasterize bars
    for ax in g.axes.flat:
        for bar in ax.patches:
            bar.set_rasterized(True)

    # Save if requested
    if save_path is not None:
        plt.tight_layout()
        g.savefig(save_path, dpi=dpi, bbox_inches="tight", format="pdf")

    plt.tight_layout()
    plt.show()


def get_cell_number_across_qc_barplots(
    doublet_file: str,
    qc_dir: str,
    datasets: List[str],
    out_path:str,
):
    """
    Generates and saves bar plots showing the number of cells retained across QC stages,
    both per sample and per dataset.

    This function:
      - Aggregates quality control (QC) statistics across samples and datasets,
        including filtering steps and doublet removal.
      - Melts the QC summary DataFrame into long format suitable for plotting.
      - Generates two bar plots:
          1. Cells per sample at each QC stage
          2. Cells per dataset at each QC stage 
      - Saves the resulting plots as PDF files to the specified output path.

    Parameters:
        doublet_file (str): Path to the CSV file containing doublet detection results.
        qc_dir (str): Directory containing per-dataset QC summary files.
        datasets (list of str): List of dataset identifiers to include in the analysis.
        out_path (str): Path to the output directory where the plots will be saved.

    Returns:
        None
    """


    qc = summarize_qc_counts(
        doublet_file = doublet_file,
        qc_dir = qc_dir,
        datasets = datasets)

    qc_melt = melt_qc_df(qc)

    plot_per_sample_qc_barplot(
        qc_melt,
        save_path=out_path + "/cells_per_sample_across_pipeline.pdf",)

    plot_per_dataset_qc_barplot(
        qc_melt,
        save_path=out_path + "/cells_per_dataset_across_pipeline.pdf",)

    

