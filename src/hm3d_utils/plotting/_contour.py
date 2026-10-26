import warnings

from os import PathLike
from typing import Literal, Mapping, Optional, Sequence

import matplotlib.axes
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import zarr

from matplotlib.ticker import ScalarFormatter, MaxNLocator
from scipy.ndimage import gaussian_filter
from skimage import measure

def plot_clusters_collapsed(
    store: zarr.Group,
    cluster_id_list: Sequence[str],
    slice_range: tuple[int, int],
    regions: Sequence[str] = ("HY",),
    modality: Literal["rel_mass", "prob"] = "rel_mass",
    view: Literal["coronal", "sagittal"] = "coronal",
    thresholds: Sequence[float] = (0.70, 0.90),
    alphas: Sequence[float] = (0.3, 0.6),
    linestyles: Sequence[Optional[Literal["solid", "dashed", "dotted"]]] = ( "dashed", "solid", ),
    linewidths: Sequence[float] = (0.8, 0.8),
    axis0_range: Optional[tuple[float, float]] = None,
    axis1_range: Optional[tuple[float, float]] = None,
    figsize: tuple[float, float] = (8, 8),
    label_fontsize: float = 7,
    label_contour: Literal["inner", "outer"] = "inner",
    region_color: str = "darkgrey",
    region_labels: bool = True,
    region_label_fontsize: float = 7,
    title: Optional[str] = None,
    prefix: str = "All clusters",
    save_path: Optional[str | PathLike[str]] = None,
    ax: Optional[matplotlib.axes.Axes] = None,
    color_mapping: Optional[ Mapping[str, str | tuple[float, float, float] | tuple[float, float, float, float]] ] = None, ) -> None:
    """
    Plot cumulative-mass contours for one or more clusters collapsed across an
    anterior–posterior (AP) slice range, optionally overlaid with anatomical
    region boundaries.

    Each cluster volume is summed across the specified AP interval and converted
    to a cumulative-mass map. Filled contours corresponding to the requested
    cumulative-mass thresholds are drawn, allowing visualization of the cluster's
    spatial extent and core. Anatomical regions can be overlaid as contour lines.

    Parameters
    ----------
    store : zarr.Group
        Zarr store containing cluster maps and region masks.

    cluster_id_list : Sequence[str]
        Cluster IDs to plot. Alternatively, descr or ABC-like lables are accepted.

    slice_range : tuple[int, int]
        AP slice interval ``(start, end)`` used to collapse the cluster volumes.
        The interval follows standard Python slicing semantics, i.e. ``end`` is
        exclusive.

    regions : Sequence[str], default=("HY",)
        Region acronyms whose boundaries should be overlaid.

    modality : {"rel_mass", "prob"}, default="rel_mass"
        Which data modality to visualize, relatie mass or membership probabilities.

    view : {"coronal", "sagittal"}, default="coronal"
        Projection used for visualization.

    thresholds : Sequence[float], default=(0.70, 0.90)
        Cumulative-mass thresholds used to generate contours. Lower thresholds
        produce larger contours, while higher thresholds highlight the cluster
        core.

    alphas : Sequence[float], default=(0.3, 0.6)
        Transparency values corresponding to each threshold.

    linestyles : Sequence[{"solid", "dashed", "dotted"} | None]
        Line styles for contour outlines. Use ``None`` to omit the outline for
        a threshold.

    linewidths : Sequence[float]
        Line widths for contour outlines. A value of ``0`` suppresses drawing
        the outline.

    axis0_range : tuple[float, float] | None, optional
        Vertical axis limits. If ``None``, the full extent is shown.

    axis1_range : tuple[float, float] | None, optional
        Horizontal axis limits. If ``None``, the full extent is shown.

    figsize : tuple[float, float], default=(8, 8)
        Figure size in inches. Ignored when an existing axis is provided.

    label_fontsize : float, default=7
        Font size of cluster labels.

    label_contour : {"inner", "outer"}, default="inner"
        Contour used for cluster label placement.

    region_color : str, default="darkgrey"
        Color used for anatomical region contours and labels.

    region_labels : bool, default=True
        Whether to display region acronym labels.

    region_label_fontsize : float, default=7
        Font size of region labels.

    title : str | None, optional
        Figure title. If ``None``, a title is generated automatically.

    prefix : str, default="All clusters"
        Prefix used when automatically generating the title.

    save_path : str | os.PathLike | None, optional
        If provided, save the figure to this location using
        ``dpi=300`` and ``bbox_inches="tight"``.

    ax : matplotlib.axes.Axes | None, optional
        Existing Matplotlib axis to draw on. If ``None``, a new figure and axis
        are created.

    color_mapping : Mapping[str, ColorLike] | None, optional
        Mapping from cluster IDs to Matplotlib-compatible colors. If ``None``,
        colors are assigned automatically from the ``tab20`` colormap.

    Returns
    -------
    None
        Draws the requested plot, optionally saves it to disk, and displays it
        when a new figure is created.
    """

    # ------------------------------------------------------------------
    # Validation & setup
    # ------------------------------------------------------------------

    assert view in ("coronal", "sagittal")
    assert label_contour in ("inner", "outer")

    cluster_id_list = [c.split()[0] for c in cluster_id_list]

    if color_mapping is None:
        palette = plt.get_cmap("tab20")
        default_colors = [
            palette((i * 2 % 20) / 20)
            for i in range(len(cluster_id_list))
        ]

    mid_ap = (slice_range[0] + slice_range[1]) // 2

    region_slices = {
        "coronal": lambda ridx: store["region"]["region_masks"][ridx, mid_ap, :, :],
        "sagittal": lambda ridx: store["region"]["region_masks"][ridx, :, :, 100].T,
    }

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=figsize)

    # ------------------------------------------------------------------
    # Plot clusters
    # ------------------------------------------------------------------

    for i, label in enumerate(cluster_id_list):

        level = f"L{label[2]}"
        idx = list(store[level]["ids"][:]).index(label)

        # --------------------------------------------------------------
        # Determine color
        # --------------------------------------------------------------

        if color_mapping is None:
            color = default_colors[i]
        else:
            if label not in color_mapping:
                warnings.warn(
                    f"Label '{label}' not found in color_mapping. "
                    "Using grey instead."
                )
            color = color_mapping.get(label, "grey")

        # --------------------------------------------------------------
        # Collapse volume
        # --------------------------------------------------------------

        volume = store[level][modality][
            idx,
            slice_range[0]:slice_range[1],
            :,
            :
        ]

        if view == "coronal":
            collapsed = np.nansum(volume, axis=0)
        else:
            collapsed = np.nansum(volume, axis=2)

        cumsum_map = np.nan_to_num(cumulative_mass_rank_map(collapsed), nan=0)

        if cumsum_map.max() == 0:
            continue

        # --------------------------------------------------------------
        # Draw contours
        # --------------------------------------------------------------

        for threshold, alpha, linestyle, linewidth in zip(
            thresholds,
            alphas,
            linestyles,
            linewidths,
        ):

            contours = measure.find_contours(cumsum_map, threshold)

            for contour in contours:

                if view == "coronal":
                    xs = contour[:, 1]
                    ys = contour[:, 0]
                else:
                    xs = contour[:, 0] + slice_range[0]
                    ys = contour[:, 1]

                ax.fill(
                    xs,
                    ys,
                    color=color,
                    alpha=alpha,
                    linewidth=0,
                )

                if linewidth and linestyle is not None:
                    ax.plot(
                        xs,
                        ys,
                        color=color,
                        linewidth=linewidth,
                        linestyle=linestyle,
                        alpha=min(alpha * 2, 1.0),
                    )

        # --------------------------------------------------------------
        # Cluster label
        # --------------------------------------------------------------

        threshold = (
            thresholds[-1]
            if label_contour == "inner"
            else thresholds[0]
        )

        contours = measure.find_contours(cumsum_map, threshold)

        if contours:

            largest = max(contours, key=len)

            if view == "coronal":
                cx = np.mean(largest[:, 1])
                cy = np.mean(largest[:, 0])
            else:
                cx = np.mean(largest[:, 0]) + slice_range[0]
                cy = np.mean(largest[:, 1])

            ax.text(
                cx,
                cy,
                label,
                fontsize=label_fontsize,
                fontweight="bold",
                color=color,
                ha="center",
                va="center",
                bbox=dict(
                    fc="white",
                    ec="none",
                    alpha=0.5,
                    pad=1,
                ),
            )

    # ------------------------------------------------------------------
    # Plot regions
    # ------------------------------------------------------------------

    SCALE_FACTOR = 3 # x3 scale for smoother region outlines

    for region in regions:

        try:
            region_idx = list(
                store["region"]["acronyms"][:]
            ).index(region)

        except ValueError:
            warnings.warn(
                f"Region '{region}' not found in region acronyms."
            )
            continue

        structure_mask2d = region_slices[view](region_idx)

        # Upscale mask by scale factor
        structure_mask2d = np.kron(structure_mask2d, np.ones((SCALE_FACTOR, SCALE_FACTOR)))
        structure_mask2d = structure_mask2d.astype(np.float32)

        # Apply Gaussian smoothing
        structure_mask2d = gaussian_filter(structure_mask2d, sigma=SCALE_FACTOR)

        contours = measure.find_contours(structure_mask2d, 0.5)
        for contour in contours:
            ax.plot(contour[:, 1] / SCALE_FACTOR, contour[:, 0] / SCALE_FACTOR, linewidth=0.8, c=region_color)

        if region_labels:

            ys, xs = np.where(structure_mask2d > 0)

            if len(xs):

                ax.text(
                    np.mean(xs) / SCALE_FACTOR,
                    np.mean(ys)/ SCALE_FACTOR,
                    region,
                    fontsize=region_label_fontsize,
                    color=region_color,
                    fontstyle="italic",
                    ha="center",
                    va="center",
                )

    # ------------------------------------------------------------------
    # Formatting
    # ------------------------------------------------------------------

    ax.set_aspect("equal")
    ax.axis("off")

    if axis0_range is None:
        ax.invert_yaxis()
    else:
        ax.set_ylim(axis0_range[1], axis0_range[0])

    if axis1_range is not None:
        ax.set_xlim(*axis1_range)

    ax.set_title(
        title or f"{prefix} | {view} | AP {slice_range[0]}-{slice_range[1]}"
    )

    if standalone:
        plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    if standalone:
        plt.show()


def plot_single_cluster(
    store: zarr.Group,
    cluster_id: str,
    slice_range: tuple[int, int],
    regions: Sequence[str] = ("HY",),
    modality: Literal["rel_mass", "prob"] = "rel_mass",
    view: Literal["coronal", "sagittal"] = "coronal",
    sagittal_slice: int = 100,
    mass_threshold: float = 0.50,
    fill_mode: Literal["flat", "density"] = "flat",
    axis0_range: Optional[tuple[float, float]] = None,
    axis1_range: Optional[tuple[float, float]] = None,
    figsize: tuple[float, float] = (6, 6),
    label_fontsize: float = 9,
    title_bold: bool = False,
    region_color: str = "darkgrey",
    region_labels: bool = True,
    region_labels_bold: bool = False,
    region_label_fontsize: float = 7,
    cluster_color: Optional[
        str
        | tuple[float, float, float]
        | tuple[float, float, float, float]
    ] = None,
    outline_linewidth: float = 1.0,
    region_linewidth: float = 0.8,
    vmin_frac: float = 0.05,
    title: Optional[str] = None,
    save_path: Optional[str | PathLike[str]] = None,
    ax: Optional[matplotlib.axes.Axes] = None,
    add_colorbar: bool = False,
    cbar_ax: Optional[matplotlib.axes.Axes] = None,
    colorbar_label_fontsize: float = 12,
    colorbar_tick_fontsize: float = 10,
    make_colorbar_bold: bool = False,
    colorbar_n_ticks: Optional[int] = None,

) -> None:
    """
    Plot a single cluster collapsed across an AP slice range.

    The cluster is identified from ``cluster_id``. The corresponding hierarchy
    level and index are automatically determined from the zarr store. The
    cluster volume is collapsed along the AP axis and visualized either as a
    uniform mask or as a density map.

    Parameters
    ----------
    store : zarr.Group
        Zarr store containing cluster maps and region information.

    cluster_id : str
        Cluster identifier. Descriptive labels are accepted; only the first
        whitespace-separated token is used.

    slice_range : tuple[int, int]
        AP slice interval ``(start, end)`` used for collapsing the cluster.
        The interval follows Python slicing semantics, i.e. ``end`` is
        exclusive.

    regions : Sequence[str], default=("HY",)
        Region acronyms whose boundaries should be overlaid.

    modality : {"rel_mass", "prob"}, default="rel_mass"
        Cluster modality to visualize.

    view : {"coronal", "sagittal"}, default="coronal"
        Projection used for visualization.

    sagittal_slice : int, default=100
        Lateral-medial slice index used for sagittal region masks.

    mass_threshold : float, default=0.5
        Fraction of cumulative cluster mass included in the displayed mask.

    fill_mode : {"flat", "density"}, default="flat"
        Visualization mode. ``flat`` displays a uniform cluster mask.
        ``density`` displays the collapsed cluster values within the mask.

    axis0_range : tuple[float, float] | None
        Vertical axis limits.

    axis1_range : tuple[float, float] | None
        Horizontal axis limits.

    figsize : tuple[float, float], default=(6, 6)
        Figure size in inches.

    label_fontsize : float, default=9
        Font size of the cluster title label.

    title_bold: float, default=False
        Whether to make title bold.

    region_labels : bool, default=True
        Whether to show anatomical region labels.

    region_label_fontsize : float, default=7
        Font size of region labels.

    region_labels_bold: bool, default=False
        Whether to plot region labels in bold.

    cluster_color : optional
        Matplotlib-compatible color for the cluster. If ``None``, the color defaults to the MPG green #116656.

    outline_linewidth : float, default=1.0
        Line width of the cluster boundary.

    region_linewidth : float, default=0.8
        Line width of anatomical region contours.

    vmin_frac : float, default=0.05
        Minimum color intensity fraction used in density mode.

    title : str | None
        Plot title. Automatically generated if ``None``.

    save_path : str | os.PathLike | None
        If provided, save the figure to this path.

    ax : matplotlib.axes.Axes | None
        Existing axis to draw on. If ``None``, a new figure is created.

    add_colorbar: bool = False
        Adds colorbar.

    cbar_ax: matplotlib.axes.Axes | None
        Existing axis to draw colorbar on. If ``None``, the cbar will be plotted on the one passed to x
            or the one created anew if ax=None.

    colorbar_label_fontsize: float = 12
        Fontsize for the cbar label.

    colorbar_tick_fontsize: float = 10
        Fontsize for the ticks of the cbar.

    make_colorbar_bold: bool = False.
        Makes all cbar labels bold.

    colorbar_n_ticks : int | None, default=None
        Approximate number of ticks on the colorbar. If ``None``, Matplotlib's
        default tick locator is used.

    Returns
    -------
    None
    """

    # ------------------------------------------------------------------
    # Validation and setup
    # ------------------------------------------------------------------

    if view not in ("coronal", "sagittal"):
        raise ValueError(
            f"view must be 'coronal' or 'sagittal', got {view!r}"
        )

    if fill_mode not in ("flat", "density"):
        raise ValueError(
            f"fill_mode must be 'flat' or 'density', got {fill_mode!r}"
        )

    cluster_id = cluster_id.split()[0]

    level_key = f"L{cluster_id[2]}"

    try:
        idx = list(store[level_key]["ids"][:]).index(cluster_id)
    except ValueError:
        raise ValueError(
            f"Cluster '{cluster_id}' not found in {level_key}."
        )

    mid_ap = (slice_range[0] + slice_range[1]) // 2

    if cluster_color is None:
        cluster_color = "#116656"

    region_slices = {
        "coronal": lambda ridx: store["region"]["region_masks"][
            ridx, mid_ap, :, :
        ],
        "sagittal": lambda ridx: store["region"]["region_masks"][
            ridx, :, :, sagittal_slice
        ].T,
    }

    # ------------------------------------------------------------------
    # Collapse volume
    # ------------------------------------------------------------------

    volume = store[level_key][modality][
        idx,
        slice_range[0]:slice_range[1],
        :,
        :,
    ]

    if view == "coronal":
        collapsed = np.nansum(volume, axis=0)
    else:
        collapsed = np.nansum(volume, axis=2)

    cumsum_map = np.nan_to_num(
        cumulative_mass_rank_map(collapsed),
        nan=0,
    )

    if cumsum_map.max() == 0:
        warnings.warn(
            f"No data for cluster '{cluster_id}' "
            "in the selected slice range."
        )
        return

    mask = cumsum_map >= (1 - mass_threshold)

        # ------------------------------------------------------------------
    # Create figure or use provided axis
    # ------------------------------------------------------------------

    standalone = ax is None

    if standalone:
        fig, ax = plt.subplots(figsize=figsize)

    # ------------------------------------------------------------------
    # Plot cluster
    # ------------------------------------------------------------------

    if fill_mode == "density":

        density = collapsed.astype(float).copy()
        density[~mask] = np.nan

        rgba = mcolors.to_rgba(cluster_color)

        light_color = (
            tuple(
                1 - vmin_frac * (1 - channel)
                for channel in rgba[:3]
            )
            + (1.0,)
        )

        cmap = mcolors.LinearSegmentedColormap.from_list(
            "cluster_cmap",
            [light_color, cluster_color],
        )

        if view == "sagittal":
            extent = [
                slice_range[0],
                slice_range[1],
                density.shape[1],
                0,
            ]

            im = ax.imshow(
                density.T,
                cmap=cmap,
                aspect="equal",
                extent=extent,
                origin="upper",
                interpolation="nearest",
            )

        else:
            im = ax.imshow(
                density,
                cmap=cmap,
                aspect="equal",
                origin="upper",
                interpolation="nearest",
            )

    else:
        contours = measure.find_contours(
            mask.astype(float),
            0.5,
        )

        for contour in contours:

            if view == "sagittal":
                xs = contour[:, 0] + slice_range[0]
                ys = contour[:, 1]
            else:
                xs = contour[:, 1]
                ys = contour[:, 0]

            ax.fill(
                xs,
                ys,
                color=cluster_color,
                alpha=0.6,
                linewidth=0,
            )

            ax.plot(
                xs,
                ys,
                color=cluster_color,
                linewidth=outline_linewidth,
            )

    # ------------------------------------------------------------------
    # Cluster outline for density mode
    # ------------------------------------------------------------------

    if fill_mode == "density":

        contours = measure.find_contours(
            mask.astype(float),
            0.5,
        )

        for contour in contours:

            if view == "sagittal":
                xs = contour[:, 0] + slice_range[0]
                ys = contour[:, 1]
            else:
                xs = contour[:, 1]
                ys = contour[:, 0]

            ax.plot(
                xs,
                ys,
                color=cluster_color,
                linewidth=outline_linewidth,
            )

    # ------------------------------------------------------------------
    # Region contours and labels
    # ------------------------------------------------------------------

    SCALE_FACTOR = 3 # x3 upscale to smooth region outlines

    for region in regions:

        try:
            region_idx = list(
                store["region"]["acronyms"][:]
            ).index(region)

        except ValueError:
            warnings.warn(
                f"Region '{region}' not found in region acronyms."
            )
            continue

        structure_mask2d = region_slices[view](region_idx)

        # Upscale mask by scale factor
        structure_mask2d = np.kron(structure_mask2d, np.ones((SCALE_FACTOR, SCALE_FACTOR)))
        structure_mask2d = structure_mask2d.astype(np.float32)

        # Apply Gaussian smoothing
        structure_mask2d = gaussian_filter(structure_mask2d, sigma=SCALE_FACTOR)

        contours = measure.find_contours(structure_mask2d, 0.5)
        for contour in contours:
            ax.plot(contour[:, 1] / SCALE_FACTOR, contour[:, 0] / SCALE_FACTOR, linewidth=0.8, c=region_color)


        if region_labels:

            ys, xs = np.where(structure_mask2d > 0)

            if len(xs) > 0:

                ax.text(
                    np.mean(xs) / SCALE_FACTOR,
                    np.mean(ys) / SCALE_FACTOR,
                    region,
                    fontsize=region_label_fontsize,
                    color=region_color,
                    ha="center",
                    va="center",
                    fontstyle="italic",
                    fontweight="bold" if region_labels_bold else None
                )

    # ------------------------------------------------------------------
    # Formatting
    # ------------------------------------------------------------------

    ax.set_aspect("equal")
    ax.axis("off")

    if axis0_range is not None:
        ax.set_ylim(
            axis0_range[1],
            axis0_range[0],
        )
    else:
        ax.invert_yaxis()

    if axis1_range is not None:
        ax.set_xlim(
            axis1_range[0],
            axis1_range[1],
        )

    ax.set_title(
        title
        or (
            f"{cluster_id} | "
            f"{view} | "
            f"AP {slice_range[0]}-{slice_range[1]} | "
            f"m={mass_threshold}"
        ),
        fontsize=label_fontsize,
        fontweight="bold" if title_bold else None
    )

    if fill_mode == "density" and add_colorbar:

        if cbar_ax is None:
            cbar = ax.figure.colorbar(
                im,
                ax=ax,
                fraction=0.046,
                pad=0.04,
            )
        else:
            cbar = ax.figure.colorbar(
                im,
                cax=cbar_ax,
            )

        cbar.set_label(
            "rel. m" if modality == "rel_mass" else "p",
            fontsize=colorbar_label_fontsize,
            fontweight="bold" if make_colorbar_bold else None,
        )

        cbar.ax.tick_params(
            labelsize=colorbar_tick_fontsize,
        )


        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((0, 0))  # always use scientific notation

        cbar.formatter = formatter

        if colorbar_n_ticks is not None:
            cbar.locator = MaxNLocator(nbins=colorbar_n_ticks)

        cbar.update_ticks()
        cbar.ax.yaxis.get_offset_text().set_fontsize(colorbar_tick_fontsize)
        offset = cbar.ax.yaxis.get_offset_text()
        offset.set_x(0.4*colorbar_tick_fontsize)  # increase to move further right
        if make_colorbar_bold:
            for tick in cbar.ax.get_yticklabels():
                tick.set_fontweight("bold")
                cbar.ax.yaxis.get_offset_text().set_fontweight("bold")

    # ------------------------------------------------------------------
    # Save / show
    # ------------------------------------------------------------------

    if standalone:
        plt.tight_layout()

    if save_path is not None:
        ax.figure.savefig(
            save_path,
            dpi=300,
            bbox_inches="tight",
        )

    if standalone:
        plt.show()



#################################
########### helper
#################################


def cumulative_mass_rank_map(A: np.ndarray) -> np.ndarray:
    """
    Compute mass-weighted cumulative distribution of a 3D array,
    ignoring NaN values. NaNs are preserved in the output.

    Parameters
    ----------
    A : np.ndarray
        Input array, can contain NaNs.

    Returns
    -------
    np.ndarray
        Array of same shape as A, with cumulative mass values.
        NaNs in the input remain NaN in the output.
    """
    flat = A.ravel()

    # Mask valid (non-NaN) values
    valid_mask = ~np.isnan(flat)
    mass = np.zeros_like(flat, dtype=float)
    mass[valid_mask] = flat[valid_mask] / np.nansum(flat)

    # Sort valid values in descending order
    order = np.argsort(mass[valid_mask])[::-1]
    cdf_sorted = 1 - np.cumsum(mass[valid_mask][order])

    # Scatter back to original positions
    cdf = np.full_like(flat, np.nan, dtype=float)
    cdf[valid_mask] = 0  # placeholder
    cdf[valid_mask] = cdf_sorted[np.argsort(order)]

    return cdf.reshape(A.shape)