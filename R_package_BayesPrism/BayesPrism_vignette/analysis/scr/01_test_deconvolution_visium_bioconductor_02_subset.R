# AIM ---------------------------------------------------------------------
# try to create a spatial experiment object from the count matrix and the image file
# the count data obtained from the deconvolution process

# libraries ---------------------------------------------------------------
library(Seurat)
library(BayesPrism)
library(tidyverse)
library(SeuratData)
library(AnnotationHub)
library(patchwork)
library(RColorBrewer)
library(scuttle)

# custom functions --------------------------------------------------------
# library(RColorBrewer)
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))

# read in data ------------------------------------------------------------
# read in the object
# save the small object
brain_small <- readRDS("../out/object/brain_small_bioconductor.rds")

plotVisium(brain_small, spots=T, point_size=1,zoom = T,show_axes = T)

ggspavis::plotSpots(brain_small,annotate = "log_sum",show_axes = T,point_size = 2) + 
  theme(legend.position = "top") +
  scale_color_viridis_c(option = "turbo")

# read in the output of bayesprims
bp_results <- readRDS("../out/object/bp_results_test_visium_bioconductor.rds")

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
  
  # convert it to integer
  storage.mode(data) <- "integer"
  
  return(data)
}) %>%
  setNames(unique(cell_types))

# notice that the number of genes matches the one present in the reference after clean up and gene selection
lapply(list_bulk, function(x){
  dim(x)
})

# dim(sc.dat.filtered.pc)

# is it of integer
is.integer(list_bulk_integer$Glutamatergic)

storage.mode(list_bulk_integer$Glutamatergic)

# sample routine to read in data ------------------------------------------
# extract the image data from the croppet object. in the case of a full dataest this step will be easier as we can simply read in the image as usual
# img <- readImgData(path = file.path(folder, "spatial"),sample_id = "sample01")
# in this case extract it from the subset object

img2 <- imgData(brain_small)

# extract the count data
# this one comes form the deconvolution processing
# counts <- assay(brain_small)
counts <- list_bulk_integer$Glutamatergic %>%
  DelayedArray()

# extract the rowData: construct observation & feature metadata
# notice that in this case the gene names might not be unique, maybe we can consider using the gene id before the deconvolution processing
LUT_genes <- data.frame(Symbol = rownames(counts))
# left_join(rowData(brain_small) %>% data.frame(),by = c("Symbol"))

rd2 <- S4Vectors::DataFrame(
  LUT_genes
  )

# extract the metadata
xyz2 <- cbind(colData(brain_small),spatialCoords(brain_small))[,c("in_tissue", "array_row", "array_col",
  "pxl_row_in_fullres", "pxl_col_in_fullres")]

# confirm the order of the barcode is the same
identical(rownames(xyz2),colnames(counts))

spe2 <- SpatialExperiment(
  assays = list(counts = counts),
  rowData = rd2, 
  colData = DataFrame(xyz2), 
  spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
  imgData = img2,
  sample_id = "sample01")

# test the plotting
plotVisium(spe2, spots=T, point_size=2,zoom = T,show_axes = T)

ggspavis::plotSpots(spe2,show_axes = T,point_size = 2) + 
  theme(legend.position = "top") +
  scale_color_viridis_c(option = "turbo")

# loop the process for all the cell types
# the image is the same for all the individual samples
img <- imgData(brain_small)

# loop the counts
# x <- "GABAergic"
data_list <- lapply(names(list_bulk_integer),function(x){
  # check the progress
  print(x)
  
  counts <- list_bulk_integer[[x]] %>%
    DelayedArray()
  
  # extract the rowData: construct observation & feature metadata
  # notice that in this case the gene names might not be unique, maybe we can consider using the gene id before the deconvolution processing
  LUT_genes <- data.frame(Symbol = rownames(counts))
  rd2 <- S4Vectors::DataFrame(
    LUT_genes
  )
  
  # extract the metadata
  xyz2 <- cbind(colData(brain_small),spatialCoords(brain_small))[,c("in_tissue", "array_row", "array_col",
                                                                    "pxl_row_in_fullres", "pxl_col_in_fullres")]
  
  # confirm the order of the barcode is the same
  identical(rownames(xyz2),colnames(counts))
  
  spe2 <- SpatialExperiment(
    assays = list(counts = counts),
    rowData = rd2, 
    colData = DataFrame(xyz2), 
    spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    imgData = img2,
    sample_id = "sample01")
  
  # add metadata
  is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe2)$Symbol)
  table(is_mito)
  
  # calculate per-spot QC metrics and store in colData
  spe2 <- addPerCellQC(spe2, subsets = list(mito = is_mito))
  
  # add the log sum metrics
  spe2$log_sum <- log(spe2$sum)
  
  return(spe2)
}) %>%
  setNames(names(list_bulk_integer))

# plot the data
list_plot <- pmap(list(data_list,names(data_list)), function(x,nm){
  # limit the color scale
  ggspavis::plotSpots(x,annotate = "log_sum",point_size = 2) + 
    ggtitle(nm) +
    theme(legend.position = "top") +
    scale_color_viridis_c(option = "turbo")
})

wrap_plots(list_plot)

list_plot2 <- pmap(list(data_list,names(data_list)), function(x,nm){
  # limit the color scale
  plotVisium(x, spots=T, point_size=2,zoom = T,show_axes = T,annotate = "log_sum") +
    ggtitle(nm)
})

wrap_plots(list_plot2)
