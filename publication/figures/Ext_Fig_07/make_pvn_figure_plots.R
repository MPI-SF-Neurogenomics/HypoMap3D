##########
### Load objects for figure plotting
##########


library(Seurat)
library(dplyr)
library(ggplot2)
library(pals,lib.loc = "/cephfs/scratch/bruening_scratch/lsteuernagel/R/user_lib/x86_64-pc-linux-gnu-library/4.4")

figure_outputpath = "/REDACTED-path/HypoMap3D_revision_figures/extended_figure_pvn/figure_plots/"
dir.create(figure_outputpath,recursive = TRUE)

# MFISH_Sim1_seurat
MFISH_Sim1_figures = readRDS("/REDACTED-path/HypoMap3D_figure_objects/MFISH_Sim1_figures.rds")

# Seq_Sim1_seurat
Seq_Sim1_figures = readRDS("/REDACTED-path/HypoMap3D_figure_objects/Seq_Sim1_figures.rds")

get_tab20_colors <- function(n) {
  base_cols <- tableau20(20)
  colorRampPalette(base_cols)(n)
}

##########
### sequencing Sim1 UMAP  plots
##########

pt_size = 0.3
label_size = 5.5

p1 = DimPlot(Seq_Sim1_figures,group.by = "seurat_clusters",reduction = "umap.rpca",label = TRUE,repel = TRUE,pt.size = pt_size,label.size = label_size,
             cols = get_tab20_colors(length(unique(Seq_Sim1_figures@meta.data$seurat_clusters)))) +NoLegend()+NoAxes()+ggtitle(NULL)
p2 = DimPlot(Seq_Sim1_figures,group.by = "L3_descr_filt",reduction = "umap.rpca",label = TRUE,repel = TRUE,pt.size = pt_size,label.size = label_size,
             cols = get_tab20_colors(length(unique(Seq_Sim1_figures@meta.data$L3_descr_filt)))) +NoLegend()+NoAxes()+ggtitle(NULL)
p3 = DimPlot(Seq_Sim1_figures,group.by = "L4_descr_filt",reduction = "umap.rpca",label = TRUE,repel = TRUE,pt.size = pt_size,label.size = label_size-1,
             cols = get_tab20_colors(length(unique(Seq_Sim1_figures@meta.data$L4_descr_filt)))) +NoLegend()+NoAxes()+ggtitle(NULL)
seq_umap_plot = cowplot::plot_grid(p1,p2,p3,nrow=1)

ggsave(filename = paste0(figure_outputpath,"seq_sim1_umap_clusters",".png"),
       plot = seq_umap_plot,
       dpi=600,
       width = 750,
       height = 250,
       device = png,
       units="mm"
       
)

##########
### MERFISH Sim1 UMAP  plots
##########

pt_size = 0.3
label_size = 5.5
p1 = DimPlot(MFISH_Sim1_figures,group.by = "seurat_clusters",reduction = "umap.rpca",label = TRUE,repel = TRUE,pt.size = pt_size,label.size = label_size,
             cols = get_tab20_colors(length(unique(MFISH_Sim1_figures@meta.data$seurat_clusters)))) +NoLegend()+NoAxes()+ggtitle(NULL)
p2 = DimPlot(MFISH_Sim1_figures,group.by = "L3_descr_filt",reduction = "umap.rpca",label = TRUE,repel = TRUE,pt.size = pt_size,label.size = label_size,
             cols = get_tab20_colors(length(unique(MFISH_Sim1_figures@meta.data$L3_descr_filt)))) +NoLegend()+NoAxes()+ggtitle(NULL)
p3 = DimPlot(MFISH_Sim1_figures,group.by = "L4_descr_filt",reduction = "umap.rpca",label = TRUE,repel = TRUE,pt.size = pt_size,label.size = label_size-1,
             cols = get_tab20_colors(length(unique(MFISH_Sim1_figures@meta.data$L4_descr_filt)))) +NoLegend()+NoAxes()+ggtitle(NULL)
merfish_umap_plot = cowplot::plot_grid(p1,p2,p3,nrow=1)

ggsave(filename = paste0(figure_outputpath,"merfish_sim1_umap_clusters",".png"),
       plot = merfish_umap_plot,
       dpi=600,
       width = 750,
       height = 250,
       device = png,
       units="mm"
       
)

##########
### Sankey plots between seq and hypomAp3d data
##########

source("HypoMap3D_revision_figures/extended_figure_pvn/helper_functions.R")

# Custom styling
seq_L3_sankey = plot_cluster_sankey(
  Seq_Sim1_figures,
  col1        = "seurat_clusters",
  col2        = "L3_descr",
  min_cells   = 50,
  flow_alpha  = 0.8,
  smooth      = 8,
  label_size = 3,
  label_color = "black",
  label_fill = "grey90",
  title       = "Sim1+ seq clusters vs HypoMap3D L3",
  legend_pos  = "right"
)+NoLegend()

seq_L3_sankey

seq_L4_sankey = plot_cluster_sankey(
  Seq_Sim1_figures,
  col1        = "seurat_clusters",
  col2        = "L4_descr",
  min_cells   = 50,
  flow_alpha  = 0.8,
  smooth      = 8,
  label_size = 2.5,
  label_color = "black",
  label_fill = "grey90",
  title       = "Sim1+ seq clusters vs HypoMap3D L4",
  legend_pos  = "right"
)+NoLegend()

seq_L4_sankey

## save

ggsave(filename = paste0(figure_outputpath,"seq_sim1_sankey_clusters_L3",".png"),
       plot = seq_L3_sankey,
       dpi=200,
       width = 150,
       height = 300,
       device = png,
       units="mm"
       
)

ggsave(filename = paste0(figure_outputpath,"seq_sim1_sankey_clusters_L4",".png"),
       plot = seq_L4_sankey,
       dpi=200,
       width = 150,
       height = 300,
       device = png,
       units="mm"
       
)

##########
### Sankey plots between merfish and hypomAp3d data
##########

# Custom styling
merfish_L3_sankey = plot_cluster_sankey(
  MFISH_Sim1_figures,
  col1        = "seurat_clusters",
  col2        = "L3_descr",
  min_cells   = 50,
  flow_alpha  = 0.8,
  smooth      = 8,
  label_size = 3,
  title       = "Sim1+ MERFISH clusters vs HypoMap3D L3",
  legend_pos  = "right"
)+NoLegend()

merfish_L3_sankey

merfish_L4_sankey = plot_cluster_sankey(
  MFISH_Sim1_figures,
  col1        = "seurat_clusters",
  col2        = "L4_descr",
  min_cells   = 50,
  flow_alpha  = 0.8,
  smooth      = 8,
  label_size = 2.5,
  title       = "Sim1+ MERFISH clusters vs HypoMap3D L4",
  legend_pos  = "right"
)+NoLegend()

merfish_L4_sankey

## save

ggsave(filename = paste0(figure_outputpath,"merfish_sim1_sankey_clusters_L3",".png"),
       plot = merfish_L3_sankey,
       dpi=200,
       width = 150,
       height = 300,
       device = png,
       units="mm"
       
)

ggsave(filename = paste0(figure_outputpath,"merfish_sim1_sankey_clusters_L4",".png"),
       plot = merfish_L4_sankey,
       dpi=200,
       width = 150,
       height = 300,
       device = png,
       units="mm"
       
)


##########
### Plot HypoMAp3d clusters in MERFISH data
##########

# all neuronal clusters with at least 200 mapped cells
target_clusters = names(table(MFISH_Sim1_figures@meta.data$L3_descr_filt))[table(MFISH_Sim1_figures@meta.data$L3_descr_filt) >= 200 & grepl("Gl|Gb",names(table(MFISH_Sim1_figures@meta.data$L3_descr_filt)))]
target_clusters

for(target_cluster in target_clusters){
  
  
  #target_cluster = "Gl3-2 PVH-SO Oxt"
  print(target_cluster)
  
  spatial_plot_data = MFISH_Sim1_figures@meta.data %>% dplyr::select(Cell_ID,L2_descr,L3_descr,L4_descr,seurat_clusters,x_rot,y_rot,image,bregma_level,Sample,Slide_ID)
  #expr_data = FetchData(MFISH_Sim1_figures,top_genes_cluster)
  #spatial_plot_data = cbind(spatial_plot_data,expr_data)
  spatial_plot_data$plot_col = "other"
  spatial_plot_data$plot_col[spatial_plot_data$L3_descr %in% target_cluster ] = "target"
  spatial_plot_data = spatial_plot_data %>% dplyr::arrange((plot_col))
  
  spatial_plot_data = spatial_plot_data %>%
    dplyr::filter(bregma_level %in% c("-0.5mm" , "-0.6mm" , "-0.7mm" , "-0.8mm",  "-0.9mm",  "-1.0mm" ,"-1.1mm"))
  
  # --- plot cluster -----------------------------------------------------------
  
  gene_plot = ggplot(spatial_plot_data, aes(x = x_rot, y = y_rot, color = plot_col)) +
    geom_point(size = 0.9, alpha = 1) +
    facet_wrap(~bregma_level,ncol = 7) +
    coord_equal() +                          # keeps spatial proportions correct
    scale_color_manual(values = c("other" = "#5A1178", "target" = "#D3E04E")) +  # "#3F0C56"
    #scale_color_manual(values = c("other" = "grey70", "target" = "purple")) +
    # scale_color_gradient(low="grey90",high = "darkred")+
    theme_minimal() +
    theme(
      axis.text  = element_blank(),          # remove axis numbers (optional)
      axis.title = element_blank(),          # remove x/y labels (optional)
      panel.grid = element_blank(),
      title = element_text(size=18),
      legend.position = "right",
      strip.text = element_text(size=15),
      panel.background = element_rect(fill = 'grey50')
    ) +
    guides(color = guide_legend(
      override.aes = list(size = 3)          # make legend points visible
    ))+NoLegend()+ggtitle(target_cluster)
  
  file_cluster_name = stringr::str_replace_all(make.names(target_cluster),"\\.","_")
  ggsave(filename = paste0(figure_outputpath,"merfish_pattern_",file_cluster_name,".png"),
         plot = gene_plot,
         dpi=600,
         width = 700,
         height = 150,
         device = png,
         units="mm"
         
  )
  ggsave(filename = paste0(figure_outputpath,"merfish_pattern_",file_cluster_name,".pdf"),
         plot = gene_plot,
         dpi=600,
         width = 700,
         height = 150,
         device = pdf,
         units="mm"
  )
}








