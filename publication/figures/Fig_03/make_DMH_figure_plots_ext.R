
##########
### load DMH seurat and overview
##########

figure_outputpath = "/REDACTED-path/HypoMap3D_revision_figures/main_figure_dmh/figure_plots/"

require(Seurat)
require(dplyr)
require(ggplot2)
require(RColorBrewer)

dmh_cluster_overview = data.table::fread("/REDACTED-path/HypoMap3D/revisions_data/dmh_cluster_overview.tsv",data.table = F)

DMH_seurat = readRDS("/REDACTED-path/HypoMap3D/revisions_data/dmh_seurat_object.rds")

color_vector_L3 = dmh_cluster_overview$color
names(color_vector_L3) = dmh_cluster_overview$L3_descr




