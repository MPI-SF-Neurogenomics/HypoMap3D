
##########
### load DMH seurat and overview
##########

figure_outputpath = "/cephfs/scratch/bruening_scratch/lsteuernagel/projects/HypoMap3D/HypoMap3D_revision_figures/main_figure_dmh/figure_plots/"

require(Seurat)
require(dplyr)
require(ggplot2)
require(RColorBrewer)

dmh_cluster_overview = data.table::fread("/cephfs/scratch/bruening_scratch/lsteuernagel/data/HypoMap3D/revisions_data/dmh_cluster_overview.tsv",data.table = F)

DMH_seurat = readRDS("/cephfs/scratch/bruening_scratch/lsteuernagel/data/HypoMap3D/revisions_data/dmh_seurat_object.rds")

color_vector_L3 = dmh_cluster_overview$color
names(color_vector_L3) = dmh_cluster_overview$L3_descr




