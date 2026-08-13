import pandas as pd

def get_clusters_in_region(
    regions,
    acronym,
    top_k=5,
    min_rel_mass=0.1,
    cum_mass_max=0.6,
    cum_mass_mode="including",
):
    """
    Extract clusters that have substantial relative mass in a brain region.

    A cluster is included if:

    1. The brain region has the highest relative mass for that cluster, OR
    2. The brain region is within the top_k regions AND
       relative_mass >= min_rel_mass AND
       cumulative relative mass satisfies the threshold.

    Parameters
    ----------
    regions : pd.DataFrame
        Long-format DataFrame with columns:
        ['cluster', 'acronym', 'relative_mass'].
    acronym : str
        Brain region acronym to query.
    top_k : int, default=5
        Maximum rank the region can have.
    min_rel_mass : float, default=0.1
        Minimum relative mass threshold.
    cum_mass_max : float, default=0.6
        Maximum cumulative relative mass threshold.
    cum_mass_mode : {"including", "excluding"}, default="including"
        How cumulative mass is evaluated:

        - "including": cumulative mass includes the current region
          (matches the original R implementation using cumsum()).
        - "excluding": cumulative mass is computed before the current
          region (i.e., excludes the current region's relative mass).

    Returns
    -------
    pd.DataFrame
        DataFrame with columns:
        ['cluster', 'region_rank', 'relative_mass', 'cum_mass'].
        The returned cum_mass corresponds to the selected
        cum_mass_mode.
    """

    if cum_mass_mode not in {"including", "excluding"}:
        raise ValueError(
            "cum_mass_mode must be either 'including' or 'excluding'."
        )

    # Work on a copy
    regions = regions.copy()

    # Rank regions within each cluster by decreasing relative mass
    regions = regions.sort_values(
        ["cluster", "relative_mass"],
        ascending=[True, False]
    )

    regions["rank"] = regions.groupby("cluster").cumcount() + 1

    # Cumulative mass including the current region (original R behavior)
    regions["cum_mass_including"] = (
        regions.groupby("cluster")["relative_mass"].cumsum()
    )

    # Cumulative mass before the current region
    regions["cum_mass_excluding"] = (
        regions["cum_mass_including"] - regions["relative_mass"]
    )

    # Select which cumulative mass definition to use
    if cum_mass_mode == "including":
        regions["cum_mass"] = regions["cum_mass_including"]
    else:
        regions["cum_mass"] = regions["cum_mass_excluding"]

    # Filter to the region of interest
    df_acro = regions[regions["acronym"] == acronym].copy()

    # Apply inclusion criteria
    cond1 = df_acro["rank"] == 1
    cond2 = (
        (df_acro["rank"] > 1)
        & (df_acro["rank"] <= top_k)
        & (df_acro["relative_mass"] >= min_rel_mass)
        & (df_acro["cum_mass"] < cum_mass_max)
    )

    result = (
        df_acro.loc[
            cond1 | cond2,
            ["cluster", "rank", "relative_mass", "cum_mass"],
        ]
        .rename(columns={"rank": "region_rank"})
        .reset_index(drop=True)
    )

    return result