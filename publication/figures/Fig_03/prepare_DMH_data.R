
##########
### load data
##########

figure_outputpath = "/REDACTED-path/HypoMap3D_revision_figures/main_figure_dmh/figure_plots/"

require(Seurat)
require(dplyr)
require(ggplot2)
require(RColorBrewer)

# load data 
hypomap_v09_controls_downsampled = readRDS("/REDACTED-path/HypoMap3D/v0.9_frozen/hypomap_v09_controls_downsampled500k.rds")
hypomap_v09_controls_downsampled = NormalizeData(hypomap_v09_controls_downsampled)

# final update cluster mapper
cluster_mapper_final = data.table::fread("/REDACTED-path/HypoMap3D-data/v0.9/8_annotation/20260730_descr_update.csv",data.table = F)

# cluster_annos
cluster_annos = hypomap_v09_controls_downsampled@meta.data %>% dplyr::distinct(L1_ID, L2_ID, L3_ID, L4_ID, L5_ID, L1_descr,L2_descr,L3_descr,L4_descr,L5_descr) 

# update names
cluster_annos = cluster_annos %>%
  dplyr::left_join(cluster_mapper_final %>% dplyr::select(L2_descr_old = old, L2_descr_new = new),by=c("L2_descr"="L2_descr_old")) %>%
  dplyr::left_join(cluster_mapper_final %>% dplyr::select(L3_descr_old = old, L3_descr_new = new),by=c("L3_descr"="L3_descr_old")) %>%
  dplyr::left_join(cluster_mapper_final %>% dplyr::select(L4_descr_old = old, L4_descr_new = new),by=c("L4_descr"="L4_descr_old")) %>%
  dplyr::left_join(cluster_mapper_final %>% dplyr::select(L5_descr_old = old, L5_descr_new = new),by=c("L5_descr"="L5_descr_old")) 
cluster_annos = cluster_annos %>% dplyr::select(L1_ID, L2_ID, L3_ID, L4_ID, L5_ID,L1_descr,L2_descr = L2_descr_new,L3_descr = L3_descr_new,L4_descr = L4_descr_new,L5_descr = L5_descr_new)

cluster_annos_long = hypomap_v09_controls_downsampled@meta.data %>%
  dplyr::distinct(cluster_id = L1_ID,cluster_name = L1_descr) %>%
  dplyr::bind_rows(
    hypomap_v09_controls_downsampled@meta.data %>%
      dplyr::distinct(cluster_id = L2_ID,cluster_name = L2_descr) 
  ) %>%
  dplyr::bind_rows(
    hypomap_v09_controls_downsampled@meta.data %>%
      dplyr::distinct(cluster_id = L3_ID,cluster_name = L3_descr) 
  ) %>%
  dplyr::bind_rows(
    hypomap_v09_controls_downsampled@meta.data %>%
      dplyr::distinct(cluster_id = L4_ID,cluster_name = L4_descr) 
  ) %>%
  dplyr::bind_rows(
    hypomap_v09_controls_downsampled@meta.data %>%
      dplyr::distinct(cluster_id = L5_ID,cluster_name = L5_descr) 
  )

## load markers
all_markers =  data.table::fread("/REDACTED-path/HypoMap3D_cell_type_marker_genes_all_levels.csv",data.table = F)
all_markers_updated = all_markers %>%
  dplyr::filter(padj < 0.001 & specificity >= 0.25) %>%
  dplyr::mutate(pct_diff = pctA - pctB) %>%
  dplyr::select(cluster = celltype_id, cluster_name = celltype_descr ,cluster_level = clustering, gene = gene_name, specificity, log2FC, pctA, pctB,pct_diff, padj,comparison) #%>%
# dplyr::filter(comparison %in% c())

# tree
# make tree
tree_structure = hypomap_v09_controls_downsampled@meta.data %>% dplyr::distinct(L1_ID,L2_ID,L3_ID,L4_ID,L5_ID) %>%
  dplyr::arrange(L1_ID,L2_ID,L3_ID,L4_ID)
tree_structure$L5 = as.character(tree_structure$L5)
#data.table::fwrite(tree_structure,"/REDACTED-path/HypoMap3D_v09x_treestructure.tsv",sep="\t")

# load spatial mapping results - rel masses
relMass_region_all = data.table::fread("/REDACTED-path/HypoMap3D-data/v0.9/9_objects/spatial/regionwise_relative_masses/regionwise_relative_masses_wo_fibers.csv",data.table = F)
colnames(relMass_region_all)
relMass_region_hypo = data.table::fread("/REDACTED-path/HypoMap3D-data/v0.9/9_objects/spatial/regionwise_relative_masses/regionwise_relative_masses_wo_fibers_or_non_HY_regions.csv",data.table = F)

# load region markers:
region_markers = data.table::fread("/REDACTED-path/HypoMap3D-data/v0.9/8_annotation/region_specific_markers/region_specific_markers_for_L3_annotation.csv",data.table = F)

# load region assign:
region_assignments_L3 = data.table::fread("/REDACTED-path/HypoMap3D-data/v0.9/9_objects/spatial/region_assignments/L3_region_assignments.csv",data.table = F)
region_assignments_L4 = data.table::fread("/REDACTED-path/HypoMap3D-data/v0.9/9_objects/spatial/region_assignments/L4_region_assignments.csv",data.table = F)

# 
# dot1 = Seurat::DotPlot(hypomap_v09_controls_downsampled,features = "Pomc",group.by = "L4_desc")
# dot1_data = dot1$data

target_cluster = c("Gb3-208","Gl3-25","Gl3-26","Gl3-28")
target_RelMass = relMass_region_all %>% dplyr::filter(celltype_id %in% target_cluster) %>%
  tidyr::gather(key="region",value="rel_mass",-celltype_id, -celltype_descr, -celltype_abc_like, -flag)

dmh_10pct = relMass_region_all %>% dplyr::filter(DMH >= 0.1 & grepl("Gl3|Gb3",celltype_id))

##########
### identify L3 & L4 DMH clusters
##########

# we filter for relative mass of 12.5% in DMH (the standard is 10%), but this way we can manually exclude some non-DMH clusters from the LHA as well as Ghrh neurons
dmh_clusters_L3 = region_assignments_L3 %>% dplyr::filter(region == "DMH" & relative_mass >= 0.125)  # 
dmh_clusters_L4 = region_assignments_L4 %>% dplyr::filter(region == "DMH")

# we manually add back Gb3-208
dmh_clusters_L3_add =  relMass_region_all %>% dplyr::filter(celltype_id %in% "Gb3-208") %>%
  tidyr::gather(key="region",value="relative_mass",-celltype_id, -celltype_descr, -celltype_abc_like, -flag) %>%
  dplyr::filter(region == "DMH") %>%
  dplyr::mutate(cum_mass = NA,region_rank = 3, flag = NA )
dmh_clusters_L3 = dplyr::bind_rows(dmh_clusters_L3,dmh_clusters_L3_add)

# update column names
dmh_cluster_overview = dmh_clusters_L3 %>%
  dplyr::left_join(cluster_annos[,c("L2_descr","L3_descr")] %>% dplyr::distinct(),by=c("celltype_descr"="L3_descr")) %>%
  dplyr::select(L3_ID = celltype_id , L3_descr= celltype_descr, L2_descr, region_rank, relative_mass, cum_mass, flag) %>% dplyr::arrange(desc(L3_ID))

# we manually add cluster Gb3-208 

# # add the L2 subclass
# dmh_cluster_overview  = dmh_cluster_overview %>% 
#   dplyr::left_join(tree_structure %>% dplyr::distinct(L2_ID,L4_ID) ,by=c("L4_ID"="L4_ID"))  %>%
#   dplyr::left_join(cluster_annos_long,by=c("L2_ID"="cluster_id")) %>%
#   dplyr::rename(L2_descr = cluster_name)

# generate a color palette
palette_name <- "Dark2"
base_colors <- brewer.pal(brewer.pal.info[palette_name, "maxcolors"], palette_name)

rgb_vals <- col2rgb(base_colors)
clean_colors <- base_colors[apply(rgb_vals, 2, function(x) max(x) - min(x) > 20 & max(x) > 40)]

set.seed(12)
n <- length(unique(dmh_cluster_overview$L3_descr))
colors_needed <- colorRampPalette(clean_colors)(n)  # interpolate if n > available colors
colors_needed[colors_needed == "#916B80"] = '#58386e' # manually change to a less grey color
dmh_cluster_overview$color_L3 <- setNames(sample(colors_needed), unique(dmh_cluster_overview$L3_descr))[dmh_cluster_overview$L3_descr]

color_vector_L3 = dmh_cluster_overview$color
names(color_vector_L3) = dmh_cluster_overview$L3_descr

#"#1B9E77" "#518C55" "#877A34" "#BD6812" "#CA611B" "#AE664D" "#58386e" "#7570B3" "#955BA7" "#B6479B" "#D6338F" "#D43A7A" "#AF5E5B" "#8A823C" "#66A61E" "#8AA716" "#AFA80D" "#D3AA05" "#DCA305" "#CA940D" "#B88515" "#A6761D"
#"#1B9E77" "#548B53" "#8D7830" "#C6650D" "#C46225" "#A6675A" "#886C8F" "#8068AE" "#A253A2" "#C43E96" "#E7298A" "#C04E69" "#997449" "#729928" "#7FA718" "#A6A80F" "#CCAA07" "#DFA504" "#CC950C" "#B98514" "#A6761D"
# "#1B9E77" "#4A8E59" "#7A7E3C" "#A96E1F" "#D95F02" "#C0632E" "#A6675A" "#58386e" "#7570B3" "#915EA8" "#AD4C9E" "#CA3A94" "#E7298A" "#C6486F" "#A66754" "#868638" "#66A61E" "#85A717" "#A6A80F" "#C6A908" "#E6AB02" "#D59D08" "#C6900F" "#B68316" "#A6761D"
 # '#7B7281' auf '#58386e' geändert

# For l2 colors
palette_name <- "Set1"
base_colors <- brewer.pal(brewer.pal.info[palette_name, "maxcolors"], palette_name)

rgb_vals <- col2rgb(base_colors)
clean_colors <- base_colors[apply(rgb_vals, 2, function(x) max(x) - min(x) > 20 & max(x) > 40)]
colors_needed_L2 <- colorRampPalette(clean_colors)(length(unique(dmh_cluster_overview$L2_descr)))

color_vector_L2 = colors_needed_L2
names(color_vector_L2) = unique(dmh_cluster_overview$L2_descr)

# manually defined L2 colors
color_mapping_L2 =  c(
  'Gl2-26 DMH Nkx2-4'= '#F781BF',
  'Gl2-25 DMH Hmx2'= '#B05B3A' ,
  'Gl2-11 PH-ant-LHA Otp.Bsx'= '#E8D430',
  'Gl2-36 DMH-LHA Vgll2'= '#58386e',
  'Gb2-64 TMv-PMv Tbx3'= '#CB6651',
  'Gb2-57 DMH Prdm13'= '#FFAF13',
  'Gb2-54 DMH Hmx2.Zfp536'= "#47A265",
  'Gb2-49 DMH Hmx2.Lef1'= '#4C71A4',
  'Gb2-43 DMH-LHA Gsx1'= "#E41A1C"
)


color_vector_L2 = color_mapping_L2

dmh_cluster_overview$color_L2 = color_vector_L2[match(dmh_cluster_overview$L2_descr,names(color_vector_L2))]


# "#E41A1C"
# "#4C71A4"
# "#F781BF"
# "#E8D430" 
# "#CB6651"
# "#7B7281" 
# "#B05B3A" 
# "#47A265"    
# "#F781BF" 

##########
### Highlight clusters on global UMAP
##########


hypomap_v09_controls_downsampled@meta.data$DMH_clusters_L2 = as.character(hypomap_v09_controls_downsampled@meta.data$L2_descr)
hypomap_v09_controls_downsampled@meta.data$DMH_clusters_L2[! hypomap_v09_controls_downsampled@meta.data$L3_ID %in%  dmh_cluster_overview$L3_ID & hypomap_v09_controls_downsampled@meta.data$L2_descr %in% dmh_cluster_overview$L2_descr] = "Non-DMH"

extra_col = c("grey50")
names(extra_col) =  "Non-DMH"
color_vector_L2_plus = c(color_vector_L2,extra_col)

umap_L2_DMH = DimPlot(hypomap_v09_controls_downsampled,group.by = "DMH_clusters_L2",raster.dpi = c(1536,1536),label = F,pt.size = 1.7,repel = TRUE,reduction = "umap")+
  NoLegend()+NoAxes()+scale_color_manual(values = color_vector_L2_plus,na.value = "grey90")

umap_L2_DMH


umap_L3_DMH = DimPlot(hypomap_v09_controls_downsampled,group.by = "L3_descr",raster.dpi = c(1536,1536),label = F,pt.size = 1.3,repel = TRUE,reduction = "umap")+
  NoLegend()+NoAxes()+scale_color_manual(values = color_vector_L3,na.value = "grey90")

umap_L3_DMH

### save plots
ggsave(filename = paste0(figure_outputpath,"umap_L2_DMH_global",".png"),
       plot = umap_L2_DMH,
       dpi=400,
       width = 300,
       height = 300,
       device = png,
       units="mm"
       
)
ggsave(filename = paste0(figure_outputpath,"umap_L2_DMH_global",".pdf"),
       plot = umap_L2_DMH,
       dpi=300,
       width = 300,
       height = 300,
       device = pdf,
       units="mm"
)

# save plots L3
ggsave(filename = paste0(figure_outputpath,"umap_L3_DMH_global",".png"),
       plot = umap_L3_DMH,
       dpi=400,
       width = 300,
       height = 300,
       device = png,
       units="mm"
       
)
ggsave(filename = paste0(figure_outputpath,"umap_L3_DMH_global",".pdf"),
       plot = umap_L3_DMH,
       dpi=300,
       width = 300,
       height = 300,
       device = pdf,
       units="mm"
)

##########
### DMH seurat subset
##########

# new version: include clusters from L2_ID class

DMH_seurat = subset(hypomap_v09_controls_downsampled, L2_descr %in% dmh_cluster_overview$L2_descr)
names(DMH_seurat@reductions)
DMH_seurat = RunUMAP(DMH_seurat,dims = 1:ncol(DMH_seurat@reductions$embed@cell.embeddings),reduction = "embed",reduction.name = "umap_DMH",reduction.key = "umap_DMH")

FeaturePlot(DMH_seurat,features = "Prlh",reduction = "umap_DMH",order = TRUE)

##########
### save DMH data
##########

data.table::fwrite(dmh_cluster_overview,"/REDACTED-path/HypoMap3D/revisions_data/dmh_cluster_overview.tsv",sep="\t")

saveRDS(DMH_seurat,"/REDACTED-path/HypoMap3D/revisions_data/dmh_seurat_object.rds")












