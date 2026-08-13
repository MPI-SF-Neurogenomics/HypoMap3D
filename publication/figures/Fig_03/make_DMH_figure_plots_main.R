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

color_vector_L3 = dmh_cluster_overview$color_L3
names(color_vector_L3) = dmh_cluster_overview$L3_descr
color_vector_L3id = dmh_cluster_overview$color_L3
names(color_vector_L3id) = dmh_cluster_overview$L3_ID
# extra_col = c("grey90")
# names(extra_col) =  "Non-DMH"
# color_vector_L3_plus = c(color_vector_L3,extra_col)

## load markers
all_markers =  data.table::fread("/REDACTED-path/HypoMap3D_cell_type_marker_genes_all_levels.csv",data.table = F)
all_markers_updated = all_markers %>%
  dplyr::filter(padj < 0.001 & specificity >= 0.25) %>%
  dplyr::mutate(pct_diff = pctA - pctB) %>%
  dplyr::select(cluster = celltype_id, cluster_name = celltype_descr ,cluster_level = clustering, gene = gene_name, specificity, log2FC, pctA, pctB,pct_diff, padj,comparison) #%>%

top_markers_all = all_markers_updated %>% dplyr::filter( cluster %in%  dmh_cluster_overview$L3_ID & comparison == "OneVsAll" ) %>%
  dplyr::filter(!grepl("Rik|Gm[0-9]",gene)) %>%
  dplyr::group_by(cluster) %>% dplyr::slice_max(n = 100,order_by = specificity) %>%
  dplyr::mutate(cluster = factor(cluster,levels =  unique(dmh_cluster_overview$L3_ID))) %>% dplyr::arrange(cluster)

# load dmh markers
# dmh_markers  = all_markers %>%
#   dplyr::filter(comparison == "DMH") %>%
#   dplyr::filter(padj < 0.001 & specificity >= 0.25) %>%
#   dplyr::mutate(pct_diff = pctA - pctB) %>%
#   dplyr::select(cluster = celltype_id, cluster_name = celltype_descr ,cluster_level = clustering, gene = gene_name, specificity, log2FC, pctA, pctB,pct_diff, padj,comparison) #%>%

# alternatively: use markers generated during annotations
dmh_markers = data.table::fread("/REDACTED-path/HypoMap3D-data/v0.9/8_annotation/region_specific_markers/region_specific_markers_for_L3_annotation.csv",data.table = F)  %>%
  dplyr::mutate(cluster_abc_like = cluster , cluster = stringr::str_extract(cluster_abc_like,"G[lb][0-9]+-[0-9]+")) %>%
  dplyr::rename( gene = gene_name) %>%
  dplyr::filter( cluster %in%  dmh_cluster_overview$L3_ID ) %>%
  dplyr::filter(!grepl("Rik|Gm[0-9]",gene)) %>%
  dplyr::group_by(cluster) %>% dplyr::slice_max(n = 100,order_by = specificity) %>%
  dplyr::select(cluster, cluster_abc_like, gene, specificity, log2FC, pctA, pctB, padj) %>%
  dplyr::mutate(cluster = factor(cluster,levels =  unique(dmh_cluster_overview$L3_ID))) %>% dplyr::arrange(cluster)


##########
### Plot DMH L3 clusters on UMAP
##########

# only id
# DMH_seurat@meta.data$DMH_clusters_L3 = NA
# DMH_seurat@meta.data$DMH_clusters_L3[DMH_seurat@meta.data$L3_ID %in%  dmh_cluster_overview$L3_ID] = as.character(DMH_seurat@meta.data$L3_ID[DMH_seurat@meta.data$L3_ID %in%  dmh_cluster_overview$L3_ID])
# DMH_seurat@meta.data$DMH_clusters_L3  = as.character(DMH_seurat@meta.data$DMH_clusters_L3 )

DMH_seurat@meta.data$DMH_clusters_L3 = as.character(DMH_seurat@meta.data$L3_descr)
DMH_seurat@meta.data$DMH_clusters_L3[! DMH_seurat@meta.data$L3_ID %in%  dmh_cluster_overview$L3_ID & DMH_seurat@meta.data$L2_descr %in% dmh_cluster_overview$L2_descr] = NA

DMH_seurat@meta.data$DMH_clusters_L3id = as.character(DMH_seurat@meta.data$L3_ID)
DMH_seurat@meta.data$DMH_clusters_L3id[! DMH_seurat@meta.data$L3_ID %in%  dmh_cluster_overview$L3_ID & DMH_seurat@meta.data$L2_descr %in% dmh_cluster_overview$L2_descr] = NA


umap_DMH_L3_subset = DimPlot(DMH_seurat,group.by = "DMH_clusters_L3",raster.dpi = c(2048,2048),raster=TRUE,label = TRUE,label.size = 6,pt.size = 3.5,repel = TRUE,na.value = "grey90",reduction = "umap_DMH")+
  NoLegend()+NoAxes()+scale_color_manual(values = color_vector_L3,na.value = "grey90")+ggtitle(NULL)

umap_DMH_L3_subset = DimPlot(DMH_seurat,group.by = "DMH_clusters_L3id",raster.dpi = c(2048,2048),raster=TRUE,label = TRUE,label.size = 8,pt.size = 3.5,repel = TRUE,na.value = "grey90",reduction = "umap_DMH")+
  NoLegend()+NoAxes()+scale_color_manual(values = color_vector_L3id,na.value = "grey90")+ggtitle(NULL)

umap_DMH_L3_subset

ggsave(filename = paste0(figure_outputpath,"umap_DMH_L3_subset",".png"),
       plot = umap_DMH_L3_subset,
       dpi=400,
       width = 300,
       height = 300,
       device = png,
       units="mm"
       
)
ggsave(filename = paste0(figure_outputpath,"umap_DMH_L3_subset",".pdf"),
       plot = umap_DMH_L3_subset,
       dpi=300,
       width = 300,
       height = 300,
       device = pdf,
       units="mm"
)



colnames(DMH_seurat@meta.data)
umap_DMH_L2_subset = DimPlot(DMH_seurat,group.by = "L2_descr",raster.dpi = c(2048,2048),raster=TRUE,label = TRUE,label.size = 8,pt.size = 3.5,repel = TRUE,na.value = "grey90",reduction = "umap_DMH")+
  NoLegend()+NoAxes()+scale_color_manual(values = color_vector_L2,na.value = "grey90")+ggtitle(NULL)

ggsave(filename = paste0(figure_outputpath,"umap_DMH_L2_subset",".png"),
       plot = umap_DMH_L2_subset,
       dpi=400,
       width = 300,
       height = 300,
       device = png,
       units="mm"
       
)
ggsave(filename = paste0(figure_outputpath,"umap_DMH_L2_subset",".pdf"),
       plot = umap_DMH_L2_subset,
       dpi=300,
       width = 300,
       height = 300,
       device = pdf,
       units="mm"
)




##########
### load ligand and receptor data
##########

# get gene metadata
gene_metadata = data.table::fread("/REDACTED-path/HypoMap3D/gene_metadata.csv",data.table = F)


# get more
### Creating a list of GPCR genes for looking at in the human datasets. 
gpcrs_mouse = "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name&format=tsv&query=%28G-protein%20coupled%20receptor%29%20AND%20(cc_scl_term%3ASL-0039)%20AND%20(organism_id%3A10090)%20NOT%20(cc_scl_term%3ASL-0191)%20AND%20(ft_transmem:*)" %>%   
  url() %>%   
  gzcon() %>%   
  readr::read_tsv() %>%   
  dplyr::filter(Reviewed == "reviewed")  
gpcrs_human = "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name&format=tsv&query=%28G-protein%20coupled%20receptor%29%20AND%20(cc_scl_term%3ASL-0039)%20AND%20(organism_id%3A9606)%20NOT%20(cc_scl_term%3ASL-0191)%20AND%20(ft_transmem:*)" %>%   
  url() %>%   
  gzcon() %>%   
  readr::read_tsv() %>%   
  dplyr::filter(Reviewed == "reviewed")  
gpcrs_human_vector = unlist(sapply(gpcrs_human$`Gene Names`,function(x){strsplit(x,split=" ")[[1]]})) 
gpcrs_mouse_vector = unlist(sapply(gpcrs_mouse$`Gene Names`,function(x){strsplit(x,split=" ")[[1]]}))

gpcrs_mouse_vector_curated = gpcrs_mouse_vector[! (grepl("Olfr|Adcy|Vmn1r|V1r|Mor|Kiaa|Sele|MOR|D3Bwg|Fbxw7",gpcrs_mouse_vector) | grepl("Or[0-9]|Gm[0-9]",gpcrs_mouse_vector)) ]

# load the latest gene anno
gene_annotations_hypomap3d = data.table::fread("/REDACTED-path/HypoMap3D/resources/gene_annotations.txt",data.table = F)
ligand_genes = unique(gene_annotations_hypomap3d$gene[gene_annotations_hypomap3d$type=="ligand/source"])
ligand_genes = c(ligand_genes,"Npvf","Th")
ligand_genes = ligand_genes[ligand_genes != "Vcam1"]


# use both:
receptor_genes = unique(c("Lepr",gene_annotations_hypomap3d$gene[gene_annotations_hypomap3d$type=="receptor"],gpcrs_mouse_vector_curated))
#receptor_genes = c(receptor_genes,"Gpr88","Gpr149","Gpr50")

tf_genes = unique(gene_metadata$`Gene symbol`[gene_metadata$is_TF])
#


##########
### Clusters to focus on
##########

# Gl3-259 , top markers are all secrted genes with bin UCMA, Stc1, Calca and has interesting set of receptors with all the acetylcholine receptors expressed 
# which and how many L4_ID cell types


# Gb3-170  Gpr50


#### sleected

# Gb3-173 DMH Kcnmb2 Brs3

# Gl3-80 Met

# Gb3-146, contains multiple Glp1r/Lepr +

# Gl3-83 , the known Pdyn/Grp ppulation, refeeding

#  Gl3-117 Npvf neurons:

# Gb3-209: Skor1 , Gpr50

# Gb3-196 Slc18a2 / Th --> change to -197

# Gl3-27


##########
### Make cluster cards fro these clusters
##########

source("HypoMap3D_revision_figures/main_figure_dmh/functions.R")

target_cluster_vector = c('Gb3-173','Gb3-146', 'Gb3-209', 'Gb3-197','Gl3-80','Gl3-117', 'Gl3-27', 'Gl3-83')

cluster_vector = unique(as.character(dmh_cluster_overview$L3_ID))

all_markers_updated_used = all_markers_updated
all_markers_updated_used$specificity = all_markers_updated_used$specificity * (1+all_markers_updated_used$pctA)

markers_dmh_used = dmh_markers
markers_dmh_used$specificity = markers_dmh_used$specificity * (1+markers_dmh_used$pctA)


cluster_interest = cluster_vector[1]
min_specificity_signaling = 1.5
min_specificity_markers = 0.5
n_signaling = 4
outputpath = paste0(figure_outputpath,"/clustercards/")
dir.create(outputpath)

# assemble data
all_cluster_data_list = list()
for(cluster_interest in cluster_vector){
  
  markers_subclass = all_markers_updated_used %>%
    dplyr::filter(padj < 0.0001 & specificity >= min_specificity_markers & cluster == cluster_interest & comparison == "OneVsSisters") %>%
    dplyr::select(gene,specificity) %>%
    dplyr::filter(!grepl("Gm|Rik|Rps|Rpl|BC|-os",gene)) %>%
    dplyr::slice_max(specificity,n=10) %>%
    dplyr::mutate(purpose = "L2_ID", type ="marker")
  markers_hy = all_markers_updated_used %>%
    dplyr::filter(padj < 0.0001 & specificity >= min_specificity_markers & cluster == cluster_interest & comparison == "OneVsAll") %>%
    dplyr::select(gene,specificity) %>%
    dplyr::filter(!grepl("Gm|Rik|Rps|Rpl|BC|-os",gene)) %>%
    dplyr::slice_max(specificity,n=10)  %>%
    dplyr::mutate(purpose = "HY", type ="marker")
  
  # dmh_markers = 
  markers_dmh = markers_dmh_used %>%
    dplyr::filter(padj < 0.0001 & specificity >= min_specificity_markers & cluster == cluster_interest ) %>%
    dplyr::select(gene,specificity) %>%
    dplyr::filter(!grepl("Gm|Rik|Rps|Rpl|BC|-os",gene)) %>%
    dplyr::slice_max(specificity,n=10)  %>%
    dplyr::mutate(purpose = "DMH", type ="marker")
  
  ligands_of_interest = all_markers_updated_used %>% dplyr::filter( cluster == cluster_interest & comparison == "OneVsAll" ) %>%
    dplyr::filter(gene %in% ligand_genes  & padj < 0.001 & specificity >= min_specificity_signaling) %>%
    dplyr::select(gene,specificity) %>%
    dplyr::slice_max(specificity,n=10) %>%
    dplyr::mutate( type ="ligands")
  if(grepl("Gb",cluster_interest)){
    nt_row = data.frame(gene = "GABA",specificity = 100, type = "ligands")
  }else{
    nt_row = data.frame(gene = "GLU",specificity = 100, type = "ligands")
  }
  if(nrow(ligands_of_interest)> 0){
    ligands_of_interest = bind_rows(nt_row,ligands_of_interest)
  }else{
    ligands_of_interest = nt_row
  }
  
  receptors_of_interest = all_markers_updated_used %>% dplyr::filter( cluster == cluster_interest & comparison == "OneVsAll" ) %>%
    dplyr::filter(gene %in% receptor_genes  & padj < 0.001 & specificity >= min_specificity_signaling) %>%
    dplyr::select(gene,specificity) %>%
    dplyr::slice_max(specificity,n=10) %>%
    dplyr::mutate( type ="receptors")
  
  # add all
  this_cluster = dplyr::bind_rows(markers_hy,markers_dmh) %>%
    dplyr::bind_rows(markers_subclass) %>%
    dplyr::bind_rows(ligands_of_interest) %>%
    dplyr::bind_rows(receptors_of_interest)
  this_cluster$cluster = cluster_interest
  
  all_cluster_data_list[[cluster_interest]] = this_cluster
}

all_cluster_data = do.call(rbind,all_cluster_data_list) %>% as.data.frame()

#make it into list
clustercard_input_list = list()
for(cluster_interest in cluster_vector){
  
  # TODO: decide order, handle duplicates (take 2nd best)
  # unique(as.character(cluster_annos$L3_descr[cluster_annos$L3_ID==cluster_interest]))
  
  this_cluster = all_cluster_data[all_cluster_data$cluster == cluster_interest,]
  this_cluster_markers = this_cluster %>%
    dplyr::filter(type == "marker") #%>%
    #dplyr::mutate(star_rating = specificity_to_dots(specificity))
  
  # split into per-purpose ranked candidate lists
  purposes   <- c("HY", "DMH", "L2_ID")
  candidates <- purrr::map(purposes, ~ this_cluster_markers %>% dplyr::filter(purpose == .x)) %>% purrr::set_names(purposes)
  
  # greedy selection: pick best unused gene per purpose
  # immediately resolve shared purpose label from this_cluster_markers
  selected   <- list()
  used_genes <- character(0)
  
  for (purp in purposes) {
    chosen     <- candidates[[purp]] %>% dplyr::filter(!gene %in% used_genes) %>% dplyr::slice(1)
    shared     <- this_cluster_markers %>% dplyr::filter(gene == chosen$gene) %>% dplyr::pull(purpose) %>% unique() %>% sort() %>% paste(collapse = "/")
    selected[[purp]] <- chosen %>% dplyr::mutate(purpose_label = shared)
    used_genes <- c(used_genes, chosen$gene)
  }
  
  marker_genes_final = dplyr::bind_rows(selected) %>%
    dplyr::select(gene, specificity, purpose_label) %>%
    dplyr::rename(purpose = purpose_label)
  marker_genes_final$gene = factor(marker_genes_final$gene,levels = as.character(used_genes))
  
  ligands_final = this_cluster %>%
    dplyr::filter(type == "ligands") %>%
    dplyr::slice_max(specificity,n=n_signaling )
  receptors_final = this_cluster %>%
    dplyr::filter(type == "receptors") %>%
    dplyr::slice_max(specificity,n=n_signaling )
  
  clustercard_input_list[[cluster_interest]] = list(
    name  =   dmh_cluster_overview$L3_descr[dmh_cluster_overview$L3_ID == cluster_interest],
    color = color_vector_L3id[cluster_interest],
    markers = marker_genes_final,
    ligands = ligands_final,
    receptors = receptors_final
  )
}

# get chosen markers for cards
dmh_markers_based_on_cards_top1 = unique(as.character(unlist(sapply(clustercard_input_list,FUN = function(x){return(as.character(x$markers$gene[1]))},simplify = F))))
dmh_markers_based_on_cards = unique(as.character(unlist(sapply(clustercard_input_list,FUN = function(x){return(as.character(x$markers$gene))},simplify = F))))

# make cards
for(cluster_interest in target_cluster_vector){
  cluster_data = clustercard_input_list[[cluster_interest]]
  # cluster_data
  clustercard_plot = build_clustercard(cluster_data,
                                       font_size_header  = 6,
                                       font_size_gene    = 5,
                                       font_size_purpose = 3.5,
                                       font_size_genes   = 5)
  
  
  
  #print( clustercard_plot)
  ggsave(filename = paste0(outputpath,"clustercard_",cluster_interest,".pdf"),
         plot = clustercard_plot,
         dpi=300,
         width = 100,
         height = 80+10*max( nrow( cluster_data$ligands),nrow( cluster_data$receptors)),
         device = "pdf",
         units="mm"
         
  )
}

##########
### Make dotplots based on genes from the cluster card analysis
##########

dmh_markers_based_on_cards_top1 = unique(as.character(unlist(sapply(clustercard_input_list,FUN = function(x){return(as.character(x$markers$gene[1]))},simplify = F))))
dmh_markers_based_on_cards = unique(as.character(unlist(sapply(clustercard_input_list,FUN = function(x){return(as.character(x$markers$gene))},simplify = F))))

DMH_seurat_onlyDMH = subset(DMH_seurat, L3_descr %in% dmh_cluster_overview$L3_descr)

Idents(DMH_seurat) = "L3_descr"
dotplot_cards_markers = DotPlot(DMH_seurat,features = unique(dmh_markers_based_on_cards_top1),assay = "RNA",scale = FALSE,dot.scale = 14,idents = dmh_cluster_overview$L3_descr)
dotplot_cards_markers = dotplot_cards_markers + theme(axis.text.x = element_text(size = 20, angle = 90),axis.text.y = element_text(size = 18))+xlab(NULL)+ylab(NULL)+scale_color_gradient(low = "grey90",high = "darkred")

dotplot_cards_markers

Idents(DMH_seurat) = "L3_descr"
dotplot_cards_markers_more = DotPlot(DMH_seurat,features = unique(dmh_markers_based_on_cards),assay = "RNA",scale = FALSE,dot.scale = 10,idents = dmh_cluster_overview$L3_descr)
dotplot_cards_markers_more = dotplot_cards_markers_more + theme_bw()  + theme(axis.text.x = element_text(size = 15, angle = 90),axis.text.y = element_text(size=18))+
  xlab(NULL)+ylab(NULL)+scale_color_gradient(low = "grey90",high = "darkred")

dotplot_cards_markers_more$data$id = factor(dotplot_cards_markers_more$data$id , levels = rev(dmh_cluster_overview$L3_descr))

dotplot_cards_markers_more


# save  dotplots
ggsave(filename = paste0(figure_outputpath,"dotplot_card_markers",".png"),
       plot = dotplot_cards_markers_more,
       dpi=400,
       width = 500,
       height = 300,
       device = png,
       units="mm"
       
)

ggsave(filename = paste0(figure_outputpath,"dotplot_card_markers",".pdf"),
       plot = dotplot_cards_markers_more,
       dpi=300,
       width = 500,
       height = 300,
       device = pdf,
       units="mm"
)

##########
### ligand / receptor exploration
##########

#### Ligands:
top_markers_ligands = all_markers_updated %>% dplyr::filter( cluster %in%  dmh_cluster_overview$L3_ID & comparison == "OneVsAll" ) %>%
  dplyr::filter(gene %in% ligand_genes  & padj < 0.001 & specificity >= 1)
# Get the DotPlot object
Idents(DMH_seurat) = "L3_descr"
always_include_ligands = c("Th","Npb","Npw","Npvf","Prlh") #" some highly specific or characteristic ligands that might not meet the minimum expression threshold we use to reduce the number of displayed genes but should be in there
dotplot_ligands <- DotPlot(DMH_seurat,features = unique(c(always_include_ligands,top_markers_ligands$gene)),assay = "RNA",scale = FALSE,dot.scale = 10,idents = dmh_cluster_overview$L3_descr)
dotplot_ligands_scaled <- DotPlot(DMH_seurat,features = unique(c(always_include_ligands,top_markers_ligands$gene)),assay = "RNA",scale = TRUE,dot.scale = 10,idents = dmh_cluster_overview$L3_descr)
# Reduce based on thresholds
dotplot_ligands_data_reduced <- dotplot_ligands$data %>%
  dplyr::group_by(features.plot) %>%
  dplyr::filter(
    #   sum(pct.exp > 12.5 & avg.exp > 0.25) >= 1,
    features.plot %in% always_include_ligands | sum(pct.exp > 20 & avg.exp > 0.5) >= 1,
    !features.plot %in% c("Efna3","Efnb2","Adam9","Col4a1","Nrg4","Sema6d","Fgf1","Ptn","Ntng1","Dlk1","Sema3a","Tafa1","Ubash3b","Sema3c","Sema3e","Itgav")
  ) %>%
  dplyr::filter(pct.exp >= 2.5) %>%
  dplyr::ungroup() %>%
  as.data.frame()
# ---- Cluster genes by avg.exp profile ----
dotplot_ligands_scaled_reduced = dotplot_ligands_scaled$data[dotplot_ligands_scaled$data$id %in% dotplot_ligands_data_reduced$id,]
expr_mat <- dotplot_ligands_scaled_reduced %>%
  tidyr::pivot_wider(id_cols = features.plot, names_from = id, values_from = avg.exp, values_fill = 0) %>%
  tibble::column_to_rownames("features.plot") %>%
  as.matrix()
# Compute hierarchical clustering order
hc <- hclust(dist(expr_mat), method = "average")
gene_order_ligands <- rownames(expr_mat)[hc$order]
# ---- Apply ordering to dotplot ----
dotplot_ligands$data <- dotplot_ligands_data_reduced
dotplot_ligands$data$id <- factor(dotplot_ligands$data$id, levels = unique(dmh_cluster_overview$L3_descr))
dotplot_ligands$data$features.plot <- factor(dotplot_ligands$data$features.plot, levels = gene_order_ligands)
# Plot
dotplot_ligands + theme(axis.text.x = element_text(size = 15, angle = 90))+xlab(NULL)+ylab(NULL)+scale_color_gradient(low = "grey90",high = "darkred")


######### ------------------ ##########


### receptors:
top_markers_receptors = all_markers_updated %>% dplyr::filter( cluster %in%  dmh_cluster_overview$L3_ID & comparison == "OneVsAll" ) %>%
  dplyr::filter(gene %in% receptor_genes  & padj < 0.001 & specificity >= 1)
# Get the DotPlot object
Idents(DMH_seurat) = "L3_descr"
always_include_receptors = c("Lepr")
dotplot_receptors <- DotPlot(DMH_seurat,features = unique(c(always_include_receptors,top_markers_receptors$gene)),assay = "RNA",scale = FALSE,dot.scale = 10,idents = dmh_cluster_overview$L3_descr)
dotplot_receptors_scaled <- DotPlot(DMH_seurat,features = unique(c(always_include_receptors,top_markers_receptors$gene)),assay = "RNA",scale = TRUE,dot.scale = 10,idents = dmh_cluster_overview$L3_descr)
# Reduce based on thresholds
dotplot_receptors_data_reduced <- dotplot_receptors$data %>%
  dplyr::group_by(features.plot) %>%
  dplyr::filter(
    features.plot %in% always_include_receptors | sum(pct.exp > 15 & avg.exp > 0.15) >= 1,
    !features.plot %in% c("Il1rap","Pde10a","Erbb4","Rorb","Abca1","F2r", "F2rl2","Fzd5","Itgav") # remove some less relevant genes for neuronal signaling
  ) %>%
  dplyr::filter(pct.exp >= 2.5) %>%
  dplyr::ungroup() %>%
  as.data.frame()
# ---- Cluster genes by avg.exp profile ----
dotplot_receptors_scaled_reduced = dotplot_receptors_scaled$data[dotplot_receptors_scaled$data$id %in% dotplot_receptors_data_reduced$id,]
expr_mat <- dotplot_receptors_scaled_reduced %>%
  tidyr::pivot_wider(id_cols = features.plot, names_from = id, values_from = avg.exp, values_fill = 0) %>%
  tibble::column_to_rownames("features.plot") %>%
  as.matrix()
# Compute hierarchical clustering order
hc <- hclust(dist(expr_mat), method = "average")
gene_order_receptors <- rownames(expr_mat)[hc$order]
# ---- Apply ordering to dotplot ----
dotplot_receptors$data <- dotplot_receptors_data_reduced
dotplot_receptors$data$id <- factor(dotplot_receptors$data$id, levels = unique(dmh_cluster_overview$L3_descr))
dotplot_receptors$data$features.plot <- factor(dotplot_receptors$data$features.plot, levels = gene_order_receptors)
# Plot
dotplot_receptors + theme(axis.text.x = element_text(size = 15, angle = 90))+xlab(NULL)+ylab(NULL)+scale_color_gradient(low = "grey90",high = "darkred")

##########
### plot ligand / receptor together in one dotplot
##########

library(ggnewscale)

# Extract data from both DotPlots
dotplot_ligands_data = dotplot_ligands$data %>% mutate(set = "ligands") %>% dplyr::filter(! "Itgav" %in% features.plot)
dotplot_receptors_data = dotplot_receptors$data %>% mutate(set = "receptors") %>% dplyr::filter(! "Itgav" %in% features.plot)

genes_ligands = as.character(unique(dotplot_ligands_data$features.plot))
genes_ligands = genes_ligands[genes_ligands != "Itgav"]
genes_receptors = as.character(unique(dotplot_receptors_data$features.plot))
genes_receptors = genes_receptors[genes_receptors != "Itgav"]

# Combine data
d_combined <- bind_rows(dotplot_ligands_data, dotplot_receptors_data)

# Order features so set A comes first, then set B
d_combined$features.plot <- factor(d_combined$features.plot, levels = c(genes_ligands, genes_receptors))

d_combined$id = factor(d_combined$id , levels = rev(dmh_cluster_overview$L3_descr))

dot_scale = 9
# Plot with two color scales
combined_signaling_dotplot = ggplot(d_combined, aes(x = features.plot, y = id)) +
  # ligands
  geom_point(data = . %>% filter(set == "ligands"), aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient(low = "grey90", high = "darkblue",  name = "Avg Exp\n(ligands)") +
  ggnewscale::new_scale_color() +  # <-- key: reset color scale
  # receptors
  geom_point(data = . %>% filter(set == "receptors"),aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient(low = "grey90", high = "darkred",  name = "Avg Exp\n(receptors)") +
  # Add vertical line to separate sets
  geom_vline(xintercept = length(genes_ligands) + 0.35, linetype = "dashed", color = "black") +
  # dot.scale
  scale_size(range = c(0, dot_scale), name = "% Expressed") +
  # cosmetics
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15, angle = 90))+
  xlab(NULL)+ylab(NULL)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),axis.text.y = element_text(size=18))


combined_signaling_dotplot



# save  dotplots
ggsave(filename = paste0(figure_outputpath,"combined_signaling_dotplot",".png"),
       plot = combined_signaling_dotplot,
       dpi=300,
       width = 600,
       height = 300,
       device = png,
       units="mm"
       
)

ggsave(filename = paste0(figure_outputpath,"combined_signaling_dotplot",".pdf"),
       plot = combined_signaling_dotplot,
       dpi=300,
       width = 600,
       height = 300,
       device = pdf,
       units="mm"
)



##########
### Other
##########




