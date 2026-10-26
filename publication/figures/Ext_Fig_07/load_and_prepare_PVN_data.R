

library(Seurat)
library(dplyr)
library(ggplot2)

##########
### Load and prepare sc-seq data
##########

# For the comparative plots in the paper we use the Sim1+ subset to recapitulate the figures from the original publication

# load original dataset
Seq_Sim1_seurat = readRDS("/REDACTED-path/HypoMap3D-data/PVN_data_CellReports/Seq_Sim1.rds")
# load projection results from HypoMap3D
Seq_Sim1_projected = anndataR::read_h5ad("/REDACTED-path/HypoMap3D-data/PVN_data_CellReports/projected/Seq_Sim1_projected.h5ad", to = "Seurat")
#Seq_Sim1_projected = NormalizeData(Seq_Sim1_projected)
Seq_Sim1_projected@meta.data$Cell_ID = rownames(Seq_Sim1_projected@meta.data)

# add projection results
Seq_Sim1_seurat@meta.data = Seq_Sim1_projected@meta.data
Seq_Sim1_seurat@reductions$embed = Seq_Sim1_projected@reductions$embed

## make cleaned up L3 and l4 columns for plotting (only show clusters with at least 30 mapped cells)
min_cells = 30
Seq_Sim1_seurat@meta.data$L3_descr = as.character(Seq_Sim1_seurat@meta.data$L3_descr)
Seq_Sim1_seurat@meta.data$L4_descr = as.character(Seq_Sim1_seurat@meta.data$L4_descr)
L3_clusters_filtered = Seq_Sim1_seurat@meta.data %>% dplyr::group_by(L3_descr) %>% dplyr::count() %>% dplyr::filter(n >= min_cells)
L3_clusters_filtered$L3_descr = as.character(L3_clusters_filtered$L3_descr)
L4_clusters_filtered = Seq_Sim1_seurat@meta.data %>% dplyr::group_by(L4_descr) %>% dplyr::count() %>% dplyr::filter(n >= min_cells)
L4_clusters_filtered$L4_descr = as.character(L4_clusters_filtered$L4_descr)

Seq_Sim1_seurat@meta.data$L3_descr_filt = NA
Seq_Sim1_seurat@meta.data$L3_descr_filt[Seq_Sim1_seurat@meta.data$L3_descr %in% L3_clusters_filtered$L3_descr] = Seq_Sim1_seurat@meta.data$L3_descr[Seq_Sim1_seurat@meta.data$L3_descr %in% L3_clusters_filtered$L3_descr]
Seq_Sim1_seurat@meta.data$L4_descr_filt = NA
Seq_Sim1_seurat@meta.data$L4_descr_filt[Seq_Sim1_seurat@meta.data$L4_descr %in% L4_clusters_filtered$L4_descr] = Seq_Sim1_seurat@meta.data$L4_descr[Seq_Sim1_seurat@meta.data$L4_descr %in% L4_clusters_filtered$L4_descr]


##########
### Load and prepare MERFISH data
##########


# For the comparative plots in the paper we use the Sim1+ subset to recapitulate the figures from the original publication

MFISH_Sim1_seurat = readRDS("/REDACTED-path/HypoMap3D-data/PVN_data_CellReports/MFISH_Sim1.rds")
MFISH_Sim1_projected = anndataR::read_h5ad("/REDACTED-path/HypoMap3D-data/PVN_data_CellReports/projected/MFISH_Sim1_projected.h5ad", to = "Seurat")
#MFISH_Sim1_projected = NormalizeData(MFISH_Sim1_projected)
MFISH_Sim1_projected@meta.data$Cell_ID = rownames(MFISH_Sim1_projected@meta.data)

# add
MFISH_Sim1_seurat@meta.data = MFISH_Sim1_projected@meta.data
MFISH_Sim1_seurat@reductions$embed = MFISH_Sim1_projected@reductions$embed

## make cleaned_up cols
MFISH_Sim1_seurat@meta.data$L3_descr = as.character(MFISH_Sim1_seurat@meta.data$L3_descr)
MFISH_Sim1_seurat@meta.data$L4_descr = as.character(MFISH_Sim1_seurat@meta.data$L4_descr)
L3_clusters_filtered = MFISH_Sim1_seurat@meta.data %>% dplyr::group_by(L3_descr) %>% dplyr::count() %>% dplyr::filter(n >= 50)
L3_clusters_filtered$L3_descr = as.character(L3_clusters_filtered$L3_descr)
L3_clusters_filtered = L3_clusters_filtered[!grepl("NN",L3_clusters_filtered$L3_descr),]
L4_clusters_filtered = MFISH_Sim1_seurat@meta.data %>% dplyr::group_by(L4_descr) %>% dplyr::count() %>% dplyr::filter(n >= 50)
L4_clusters_filtered$L4_descr = as.character(L4_clusters_filtered$L4_descr)
L4_clusters_filtered = L4_clusters_filtered[!grepl("NN",L4_clusters_filtered$L4_descr),]

MFISH_Sim1_seurat@meta.data$L3_descr_filt = NA
MFISH_Sim1_seurat@meta.data$L3_descr_filt[MFISH_Sim1_seurat@meta.data$L3_descr %in% L3_clusters_filtered$L3_descr] = MFISH_Sim1_seurat@meta.data$L3_descr[MFISH_Sim1_seurat@meta.data$L3_descr %in% L3_clusters_filtered$L3_descr]
MFISH_Sim1_seurat@meta.data$L4_descr_filt = NA
MFISH_Sim1_seurat@meta.data$L4_descr_filt[MFISH_Sim1_seurat@meta.data$L4_descr %in% L4_clusters_filtered$L4_descr] = MFISH_Sim1_seurat@meta.data$L4_descr[MFISH_Sim1_seurat@meta.data$L4_descr %in% L4_clusters_filtered$L4_descr]

# Get coordinates for all images
all_coords <- lapply(names(MFISH_Sim1_seurat@images), function(img) {
  coords <- GetTissueCoordinates(MFISH_Sim1_seurat,
                                 image = img,
                                 which = "centroids")
  coords$image <- img  # track which image each cell came from
  return(coords)
}) %>%
  dplyr::bind_rows()

# add to meta
MFISH_Sim1_seurat@meta.data$Cell_ID = rownames(MFISH_Sim1_seurat@meta.data)
temp_meta = MFISH_Sim1_seurat@meta.data %>% dplyr::left_join(all_coords,by=c("Cell_ID"="cell")) %>% as.data.frame()
rownames(temp_meta) = temp_meta$Cell_ID
MFISH_Sim1_seurat@meta.data = temp_meta

##########
### Make rotations of MERFIS data
##########

# Plot each slice separately to assess rotation needed
lapply(unique(MFISH_Sim1_seurat$image), function(b) {

  df <- MFISH_Sim1_seurat@meta.data %>% filter(image == b)

  print(ggplot(df, aes(x = x, y = y, color = seurat_clusters)) +
          geom_point(size = 0.3) +
          coord_equal() +
          ggtitle(paste("Image:", b)) +
          theme_minimal()+
          NoLegend())
})

# Named vector of rotation angles per image (in degrees)
# manually assesed necessary rotation
rotation_angles <- c(
  "PVH1Region3z3L" = -110,  ## REMOVE?
  "PVH2Region1z3L" = 45,
  "PVH2Region2z3L" = 45,
  "PVH2Region3z3L" = 45,  ## REMOVE?
  "PVH2Region4z3L" = 45,
  "PVH2Region5z3L" = -131,
  "PVH2Region6z3L" = -131,
  "PVH2Region7z3L" = 45,
  "PVH4Region0z3L" = -12,
  "PVH4Region1z3L" = 165,
  "PVH4Region2z3L" = 165,
  "PVH4Region3z3L" = 165,
  "PVH4Region4z3L" = 165,
  "PVH4Region5z3L" = -12,
  "PVH4Region6z3L" = -12,
  "PVH4Region7z3L" = 170,
  "PVH4Region8z3L" = 173,
  "PVH4Region9z3L" = 170,  ## REMOVE?
  "PVH5Region0z3L" = -7,  ## REMOVE?
  "PVH5Region1z3L" = -12,  ## REMOVE?
  "PVH5Region2z3L" = 175,
  "PVH5Region3z3L" = 173,
  "PVH5Region4z3L" = 170,
  "PVH5Region5z3L" = 175,
  "PVH5Region6z3L" = 175,  ## REMOVE?
  "PVH5Region7z3L" = 175,
  "PVH5Region8z3L" = 175,
  "PVH6Region1z3L" = -158,
  "PVH6Region2z3L" = -158,
  "PVH6Region3z3L" = -158,
  "PVH6Region4z3L" = -164,
  "PVH6Region5z3L" = 20,
  "PVH6Region6z3L" = 25,
  "PVH6Region7z3L" = -158,
  "PVH6Region9z3L" = -158,
  "PVH7Region0z3L" = 108,
  "PVH7Region1z3L" = 108,
  "PVH7Region2z3L" = 108,
  "PVH7Region3z3L" = -75,
  "PVH7Region4z3L" = -76,
  "PVH7Region5z3L" = -79
)
#

# Rotation function
rotate_coords <- function(x, y, angle_deg) {
  angle_rad <- angle_deg * pi / 180

  x_rot <- x * cos(angle_rad) - y * sin(angle_rad)
  y_rot <- x * sin(angle_rad) + y * cos(angle_rad)

  return(data.frame(x_rot = x_rot, y_rot = y_rot))
}

# Apply rotation per image
metadata <- MFISH_Sim1_seurat@meta.data

metadata <- metadata %>%
  dplyr::group_by(image) %>%
  dplyr::mutate(
    # Center coordinates first before rotating
    x_centered = x - mean(x),
    y_centered = y - mean(y),
    # Apply rotation using the angle for this image
    x_rot = rotate_coords(x_centered, y_centered, rotation_angles[unique(image)])$x_rot,
    y_rot = rotate_coords(x_centered, y_centered, rotation_angles[unique(image)])$y_rot
  ) %>%
  ungroup() %>% as.data.frame()

# Add back to seurat object
MFISH_Sim1_seurat <- AddMetaData(MFISH_Sim1_seurat, metadata = metadata[, c("x_rot", "y_rot")])

# Quick check
p1 = ggplot(MFISH_Sim1_seurat@meta.data[MFISH_Sim1_seurat@meta.data$image == "PVH7Region2z3L",], aes(x = x_rot, y = y_rot, color = seurat_clusters)) +
  geom_point(size = 0.3) +
  coord_equal() +
  theme_minimal()+
  NoLegend()

p2 = ggplot(MFISH_Sim1_seurat@meta.data[MFISH_Sim1_seurat@meta.data$image == "PVH7Region2z3L",], aes(x = x, y = y, color = seurat_clusters)) +
  geom_point(size = 0.3) +
  coord_equal() +
  theme_minimal()+
  NoLegend()

cowplot::plot_grid(p1,p2,nrow = 2)

#check all
lapply(unique(MFISH_Sim1_seurat$image), function(b) {

  df <- MFISH_Sim1_seurat@meta.data %>% filter(image == b)

  print(ggplot(df, aes(x = x_rot, y = y_rot, color = seurat_clusters)) +
          geom_point(size = 0.3) +
          coord_equal() +
          ggtitle(paste("Image:", b)) +
          theme_minimal()+
          NoLegend())
})


##########
###Save objects for figure plotting
##########


# save result

# MFISH_Sim1_seurat
saveRDS(MFISH_Sim1_seurat,"/REDACTED-path/HypoMap3D-data/PVN_data_CellReports/hypomap3d_figure_objects/MFISH_Sim1_figures.rds")

# Seq_Sim1_seurat
saveRDS(Seq_Sim1_seurat,"/REDACTED-path/HypoMap3D-data/PVN_data_CellReports/hypomap3d_figure_objects/Seq_Sim1_figures.rds")
