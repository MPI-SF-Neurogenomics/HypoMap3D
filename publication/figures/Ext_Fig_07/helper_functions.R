
##########
### sankey plots
##########

# ggsankey based Sankey plot function working directly on seurat object metadata columns (generated with help of Claude)

library(ggsankey,lib.loc = "/cephfs/scratch/bruening_scratch/lsteuernagel/R/user_lib/x86_64-pc-linux-gnu-library/4.4")
library(ggplot2)
library(dplyr)
library(RColorBrewer)

#' Plot a Sankey diagram comparing two clusterings from a Seurat object
#'
#' @param seurat_obj     A Seurat object
#' @param col1           Name of the first clustering column in meta.data
#' @param col2           Name of the second clustering column in meta.data
#' @param min_cells      Minimum number of cells required to draw a connection (default: 20)
#' @param flow_alpha     Transparency of the flow ribbons (default: 0.6)
#' @param node_color     Outline colour of the nodes (default: "grey30")
#' @param smooth         Smoothness of the flow curves, higher = smoother (default: 6)
#' @param label_size     Size of the node labels (default: 3.5)
#' @param label_color    Colour of the node label text (default: "white")
#' @param label_fill     Fill colour behind the node labels (default: "grey30")
#' @param palette        A ColorBrewer palette name or a named/unnamed character vector of custom colours (default: "Set2")
#' @param title          Plot title (default: "Cluster Comparison")
#' @param base_size      Base font size for theme_minimal (default: 13)
#' @param legend_pos     Legend position, e.g. "bottom", "right", "none" (default: "bottom")
#'
#' @return A ggplot object
#'
plot_cluster_sankey <- function(
    seurat_obj,
    col1,
    col2,
    min_cells   = 20,
    flow_alpha  = 0.6,
    node_color  = "grey30",
    smooth      = 6,
    label_size  = 3.5,
    label_color = "white",
    label_fill  = "grey30",
    palette     = "Set2",
    title       = "Cluster Comparison",
    base_size   = 13,
    legend_pos  = "bottom"
) {
  
  # --- Input validation ---
  stopifnot(
    "seurat_obj must be a Seurat object"    = inherits(seurat_obj, "Seurat"),
    "col1 must be in seurat_obj@meta.data"  = col1 %in% colnames(seurat_obj@meta.data),
    "col2 must be in seurat_obj@meta.data"  = col2 %in% colnames(seurat_obj@meta.data),
    "min_cells must be a positive integer"  = is.numeric(min_cells) && min_cells >= 1,
    "flow_alpha must be between 0 and 1"    = is.numeric(flow_alpha) && flow_alpha >= 0 && flow_alpha <= 1
  )
  
  # --- Summarise cluster overlaps ---
  cluster_summary <- seurat_obj@meta.data |>
    dplyr::select(cluster1 = dplyr::all_of(col1),
                  cluster2 = dplyr::all_of(col2)) |>
    dplyr::mutate(
      cluster1 = as.character(cluster1),
      cluster2 = as.character(cluster2)
    ) |>
    dplyr::count(cluster1, cluster2, name = "n_cells") |>
    dplyr::filter(n_cells >= min_cells)
  
  if (nrow(cluster_summary) == 0) {
    stop("No connections remain after applying min_cells = ", min_cells,
         ". Try lowering the threshold.")
  }
  
  n_nodes_removed <- seurat_obj@meta.data |>
    dplyr::select(cluster1 = dplyr::all_of(col1),
                  cluster2 = dplyr::all_of(col2)) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) |>
    dplyr::count(cluster1, cluster2) |>
    dplyr::filter(n < min_cells) |>
    nrow()
  
  if (n_nodes_removed > 0) {
    message(n_nodes_removed, " connection(s) removed due to min_cells = ", min_cells)
  }
  
  # --- Optimise node order to minimise crossing connections ---
  # First pass: order col2 nodes by weighted mean position of their col1 connections
  col2_order <- cluster_summary |>
    dplyr::group_by(cluster2) |>
    dplyr::summarise(
      mean_c1 = stats::weighted.mean(as.numeric(factor(cluster1)), n_cells),
      .groups = "drop"
    ) |>
    dplyr::arrange(mean_c1) |>
    dplyr::pull(cluster2)
  
  # Second pass: order col1 nodes by weighted mean position of their col2 connections
  col1_order <- cluster_summary |>
    dplyr::mutate(cluster2_rank = as.numeric(factor(cluster2, levels = col2_order))) |>
    dplyr::group_by(cluster1) |>
    dplyr::summarise(
      mean_c2 = stats::weighted.mean(cluster2_rank, n_cells),
      .groups = "drop"
    ) |>
    dplyr::arrange(mean_c2) |>
    dplyr::pull(cluster1)
  
  # --- Build a consistent discrete colour palette across all nodes ---
  all_nodes <- union(col1_order, col2_order)
  n_nodes   <- length(all_nodes)
  
  if (length(palette) == 1 && is.character(palette)) {
    max_brewer <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
    if (n_nodes <= max_brewer) {
      node_colours <- RColorBrewer::brewer.pal(max(3, n_nodes), palette)[seq_len(n_nodes)]
    } else {
      node_colours <- grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(max_brewer, palette))(n_nodes)
    }
  } else {
    node_colours <- rep_len(palette, n_nodes)
  }
  names(node_colours) <- all_nodes
  
  # --- Reshape to ggsankey long format ---
  sankey_data <- cluster_summary |>
    dplyr::mutate(
      cluster1 = as.character(cluster1),
      cluster2 = as.character(cluster2)
    ) |>
    ggsankey::make_long(cluster1, cluster2, value = n_cells)
  
  # --- Explicitly enforce node order on the sankey_data after make_long ---
  # ggsankey uses the factor levels of `node` within each x position to
  # determine vertical order — this is the key step that was missing before
  sankey_data <- sankey_data |>
    dplyr::mutate(
      node = dplyr::case_when(
        x == "cluster1" ~ factor(node, levels = rev(col1_order)),
        x == "cluster2" ~ factor(node, levels = rev(col2_order))
      ),
      next_node = dplyr::case_when(
        next_x == "cluster1" ~ factor(next_node, levels = rev(col1_order)),
        next_x == "cluster2" ~ factor(next_node, levels = rev(col2_order)),
        is.na(next_x)        ~ factor(next_node)
      )
    )
  
  # --- Plot ---
  p <- ggplot2::ggplot(
    sankey_data,
    ggplot2::aes(
      x         = x,
      next_x    = next_x,
      node      = node,
      next_node = next_node,
      fill      = node,
      value     = value,
      label     = node
    )
  ) +
    ggsankey::geom_sankey(
      flow.alpha = flow_alpha,
      node.color = node_color,
      smooth     = smooth
    ) +
    ggsankey::geom_sankey_label(
      size  = label_size,
      color = label_color,
      fill  = label_fill
    ) +
    ggplot2::scale_fill_manual(values = node_colours) +
    ggplot2::scale_x_discrete(
      labels = c(col1, col2),
      expand = c(0.1, 0.1)
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = paste0("Connections with < ", min_cells, " cells removed"),
      x        = NULL,
      fill     = "Cluster / Annotation"
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position = legend_pos,
      axis.text.y     = ggplot2::element_blank(),
      panel.grid      = ggplot2::element_blank()
    )
  
  return(p)
}