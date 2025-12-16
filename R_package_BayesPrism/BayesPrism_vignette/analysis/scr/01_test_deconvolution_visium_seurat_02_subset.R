# AIM ---------------------------------------------------------------------
# try to create a seurat object from the count matrix and the image file
# the count data obtained from the deconvolution process

# libraries ---------------------------------------------------------------
library(Seurat)
library(BayesPrism)
library(tidyverse)
library(SeuratData)
library(AnnotationHub)
library(patchwork)
library(RColorBrewer)

# custom functions --------------------------------------------------------
# library(RColorBrewer)
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))

# read in data ------------------------------------------------------------
# import the sample reference dataset
seurat_ref <- readRDS("../data/allen_cortex.rds")

# import  the full spatial dataset
brain <- LoadData("stxBrain", type = "anterior1")
# plot to check the object
SpatialFeaturePlot(brain, features = "nCount_Spatial")

# regenerate the smaller dataset uset for the test
brain <- AddMetaData(object = brain,metadata = GetTissueCoordinates(brain))

# Define the boundaries for a rectangular region
min_x <- 5500  # Minimum imagecol value
max_x <- 6500 # Maximum imagecol value
min_y <- 3500  # Minimum imagerow value
max_y <- 5500 # Maximum imagerow value

# subset the image
brain_small <- subset(brain, subset = 
                        # Filter by X coordinate (imagecol)
                        y >= min_y & y <= max_y &
                        x >= min_x & x <= max_x)

SpatialFeaturePlot(brain_small, features = "nCount_Spatial")


# read in the output of bayesprims
bp_results <- readRDS("../out/object/bp_results_test_visium.rds")

# wrangling ---------------------------------------------------------------

# get the fraction estimated
theta <- get.fraction(bp=bp_results,
                      which.theta="final",
                      state.or.type="type")
theta

# get the expression data per celltype
# x <- "Endothelial"

# pull the cell types from the object
cell_types <- theta %>%
  colnames()

# extract the reads data per cell type after cell type deconvolution
list_bulk <- lapply(unique(cell_types), function(x){
  mat <- get.exp(bp=bp_results,
                 state.or.type="type",
                 cell.name=x)
  
  data <- mat %>%
    t()
  
  return(data)
}) %>%
  setNames(unique(cell_types))

# generate integer estimates
list_bulk_integer <- lapply(unique(cell_types), function(x){
  mat <- get.exp(bp=bp_results,
                 state.or.type="type",
                 cell.name=x)
  
  # generate integer level estimates
  # https://github.com/Danko-Lab/BayesPrism/issues/121
  # round up the Z matrix to the nearest integer and use it as input for DESeq2 (without any normalization step).
  data <- mat %>%
    t() %>%
    round()
  return(data)
}) %>%
  setNames(unique(cell_types))

# notice that the number of genes matches the one present in the reference after clean up and gene selection
lapply(list_bulk, function(x){
  dim(x)
})

# dim(sc.dat.filtered.pc)

# sample routine to read in data ------------------------------------------
# extract the image data from the croppet object. in the case of a full dataest this step will be easier as we can simply read in the image as usula
# img <- Read10X_Image(image.dir = "../data/GSE210616_RAW/", image.name = "tissue_hires_image.png")
# in this case extract it from the subset object

# single sample
img <- brain_small@images$anterior1

# confirm the assay is Spatial
Seurat::DefaultAssay(object = img)

# read in count data from the individual cell types and convert it to sparse matrix
counts <- list_bulk_integer$Glutamatergic %>%
  as.sparse()

# create the object
data <- Seurat::CreateSeuratObject(
  counts = counts , 
  project = 'test', 
  assay = 'Spatial')

# add metadata
data$slice <- 1 
data$region <- 'test' 

# confirm the barcodes are matching
# colnames(x = data)
# img@boundaries$centroids@cells
identical(colnames(x = data),img@boundaries$centroids@cells)

# add the image to the counts
data[['image']] <- img  

# sample plotting
SpatialFeaturePlot(data, features = "nCount_Spatial")

# loop the process for all the cell types
# the image is the same for all the individual samples
img <- brain_small@images$anterior1

# loop the counts
# x <- "GABAergic"
data_list <- lapply(names(list_bulk_integer),function(x){
  # check the progress
  print(x)
  
  # confirm the assay is Spatial
  Seurat::DefaultAssay(object = img)
  
  # read in count data from the individual cell types and convert it to sparse matrix
  counts <- list_bulk_integer[[x]] %>%
    as.sparse()
  
  # create the object
  data <- Seurat::CreateSeuratObject(
    counts = counts , 
    project = 'test', 
    assay = 'Spatial')
  
  # add metadata
  data$slice <- 1 
  data$region <- 'test' 
  
  # add the image to the counts
  data[['image']] <- img  
  
  return(data)
}) %>%
  setNames(names(list_bulk_integer))

# plot the data
list_plot <- pmap(list(data_list,names(data_list)), function(x,nm){
  # limit the color scale
  limit_min <- 0
  limit_max <- 7000
  round(limit_max - (limit_max-limit_min)/2)
  
  SpatialFeaturePlot(x, features = "nCount_Spatial") +
    ggtitle(nm) +
    scale_fill_gradientn(colours = SpatialColors(n = 100),
                         breaks = c(limit_min,
                                    round(limit_max - (limit_max-limit_min)/2),
                                    limit_max),
                         # labels = c("Low", "Medium", "High"),
                         limits = c(limit_min, limit_max),
                         # this one allows to cap the value above a certain threshold
                         oob = scales::squish)+
    theme(legend.position = "right")
})

wrap_plots(list_plot)
