library(tidyverse)
library(patchwork)

# ─────────────────────────────────────────────
# Plot 1: Header bar with cluster name
# ─────────────────────────────────────────────
plot_header <- function(cluster_data, font_size_header = 4.5) {
  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = cluster_data$color) +
    annotate("text", x = 0.5, y = 0.5, label = cluster_data$name, color = "white", fontface = "bold", size = font_size_header) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

# ─────────────────────────────────────────────
# Plot 2: Marker genes with specificity bars
# Bar fills from 0 to min(specificity, maxSpec)
# Gene name on left, purpose label on right
# ─────────────────────────────────────────────
plot_markers <- function(cluster_data, font_size_gene = 3, font_size_purpose = 2.5, bar_height = 0.4, maxSpec = 8) {
  
  df <- cluster_data$markers %>%
    dplyr::mutate(
      spec_capped = pmin(specificity, maxSpec),
      gene        = factor(gene, levels = rev(gene))
    )
  
  # x layout: gene label | bar (bar_start to bar_end) | purpose label
  bar_start <- 0
  bar_end   <- 1
  label_x   <- bar_start - 0.05
  purpose_x <- bar_end   + 0.05
  
  ggplot(df) +
    # background bar (full width = maxSpec)
    geom_rect(aes(xmin = bar_start, xmax = bar_end, ymin = as.numeric(gene) - bar_height/2, ymax = as.numeric(gene) + bar_height/2), fill = "gray90") +
    # filled bar (scaled to specificity)
    geom_rect(aes(xmin = bar_start, xmax = bar_start + (spec_capped / maxSpec) * (bar_end - bar_start), ymin = as.numeric(gene) - bar_height/2, ymax = as.numeric(gene) + bar_height/2), fill = cluster_data$color) +
    # gene name
    geom_text(aes(x = label_x, y = as.numeric(gene), label = gene), hjust = 1, fontface = "bold", size = font_size_gene) +
    # purpose label
    geom_text(aes(x = purpose_x, y = as.numeric(gene), label = purpose), hjust = 0, color = "gray50", size = font_size_purpose) +
    scale_x_continuous(limits = c(-1.2, 2)) +
    scale_y_continuous(limits = c(0.5, nrow(df) + 0.5), expand = c(0, 0)) +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(4, 4, 4, 4))
}

# ─────────────────────────────────────────────
# Plot 3: Ligand and receptor boxes
# Both columns start from top and grow downward
# ─────────────────────────────────────────────
plot_ligand_receptor <- function(cluster_data, font_size_genes = 3, y_step = 0.18, n_signaling = 4) {
  
  ligands   <- cluster_data$ligands   %>% dplyr::arrange(dplyr::desc(specificity)) %>% dplyr::pull(gene)
  receptors <- cluster_data$receptors %>% dplyr::arrange(dplyr::desc(specificity)) %>% dplyr::pull(gene)
  
  # pad to n_signaling so all cards have same height
  ligands   <- c(ligands,   rep("", max(0, n_signaling - length(ligands))))
  receptors <- c(receptors, rep("", max(0, n_signaling - length(receptors))))
  
  n_max  <- n_signaling
  y_top  <- 0.95
  y_bot  <- y_top - (n_max * y_step) - 0.1
  
  df_lig <- data.frame(x = 0.25, y = seq(y_top - 0.1, by = -y_step, length.out = length(ligands)),   label = ligands)
  df_rec <- data.frame(x = 0.75, y = seq(y_top - 0.1, by = -y_step, length.out = length(receptors)), label = receptors)
  
  ggplot() +
    annotate("rect", xmin = 0.02, xmax = 0.48, ymin = y_bot, ymax = y_top, fill = scales::alpha(cluster_data$color, 0.08), color = cluster_data$color, linewidth = 0.5) +
    annotate("rect", xmin = 0.52, xmax = 0.98, ymin = y_bot, ymax = y_top, fill = scales::alpha(cluster_data$color, 0.08), color = cluster_data$color, linewidth = 0.5) +
    geom_text(data = df_lig, aes(x = x, y = y, label = label), size = font_size_genes, fontface = "bold", hjust = 0.5) +
    geom_text(data = df_rec, aes(x = x, y = y, label = label), size = font_size_genes, fontface = "bold", hjust = 0.5) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(y_bot - 0.02, y_top + 0.02), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(2, 2, 2, 2))
}

# ─────────────────────────────────────────────
# Assemble: stack all three plots using patchwork
# ─────────────────────────────────────────────
build_clustercard <- function(cluster_data,
                              font_size_header  = 4.5,
                              font_size_gene    = 3,
                              font_size_purpose = 2.5,
                              font_size_genes   = 3,
                              bar_height        = 0.4,
                              maxSpec           = 8,
                              y_step            = 0.18,
                              n_signaling       = 4,
                              heights           = NULL) {
  
  p_header  <- plot_header(cluster_data, font_size_header = font_size_header)
  p_markers <- plot_markers(cluster_data, font_size_gene = font_size_gene, font_size_purpose = font_size_purpose, bar_height = bar_height, maxSpec = maxSpec)
  p_boxes   <- plot_ligand_receptor(cluster_data, font_size_genes = font_size_genes, y_step = y_step, n_signaling = n_signaling)
  
  if (is.null(heights)) {
    heights <- c(0.5, nrow(cluster_data$markers) * 0.3, n_signaling * 0.25)
  }
  
  assembled <- p_header / p_markers / p_boxes +
    plot_layout(heights = unit(heights, "null")) &
    theme(plot.background = element_blank(), plot.margin = margin(0, 0, 0, 0))
  
  assembled +
    plot_annotation(theme = theme(
      plot.background = element_rect(color = cluster_data$color, linewidth = 2, fill = "white"),
      plot.margin     = margin(0, 0, 0, 0)
    ))
}

# # ─────────────────────────────────────────────
# # Example input and usage
# # ─────────────────────────────────────────────
# cluster_data <- list(
#   name  = "PV Interneuron",
#   color = "#E05C5C",
#   markers = data.frame(
#     gene        = c("PVALB", "SST", "KCNC1"),
#     specificity = c(8.2, 2.1, 5.6),
#     purpose     = c("vs Subclass", "vs Region", "vs Nucleus")
#   ),
#   ligands = data.frame(
#     gene        = c("GABA", "GLU"),
#     specificity = c(9.1, 3.2)
#   ),
#   receptors = data.frame(
#     gene        = c("GABRA1", "KCNC1", "GRM1"),
#     specificity = c(7.2, 5.1, 2.3)
#   )
# )
# 
# build_clustercard(cluster_data)