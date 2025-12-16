# AIM  --------------------------------------------------------------------
# sample build a SpatialExperiment from scratch

# libraries ---------------------------------------------------------------
library(Seurat)
library(BayesPrism)
library(tidyverse)
library(SeuratData)
library(AnnotationHub)
library(VisiumIO)
library(ggspavis)
library(DropletUtils)

# read in the data --------------------------------------------------------
# folder <- "../data/sample_Visium/results/spaceranger/merged/V1_Mouse_Brain_Sagittal_Anterior_Section_1/"
# brain <- TENxVisium(
#   spacerangerOut = folder,
#   processing = "filtered", 
#   format="h5", 
#   images="hires") %>%
#   import()
# 
# # save the sample object
# base::saveRDS(brain,"../out/object/brain_sample_bioconductor.rds")

brain <- readRDS("../out/object/brain_sample_bioconductor.rds")

# confirm the object by producign a plots
plotVisium(brain, spots=T, point_size=1)

# explore the object
getImg(brain)
spatialCoords(brain)
imgData(brain)

# manually construct the object -------------------------------------------
# locate the raw folder
folder <- "../data/sample_Visium/results/spaceranger/merged/V1_Mouse_Brain_Sagittal_Anterior_Section_1/outs"

# read in counts directly from the output folder
fnm <- file.path(folder, "filtered_feature_bc_matrix/")
sce <- DropletUtils::read10xCounts(fnm)

# read in image data
img <- readImgData(
  path = file.path(folder, "spatial"),
  sample_id = "sample01")

# read in spatial coordinates
fnm <- file.path(folder, "spatial", "tissue_positions.csv")
xyz <- read.csv(fnm,
                header = T,
                col.names = c(
                  "barcode", "in_tissue", "array_row", "array_col",
                  "pxl_row_in_fullres", "pxl_col_in_fullres")) %>%
  # filter only the in tissue data
  filter(in_tissue == 1)

# construct observation & feature metadata
rd <- S4Vectors::DataFrame(
  symbol = rowData(sce)$Symbol)

# construct 'SpatialExperiment'
(spe <- SpatialExperiment(
  assays = list(counts = assay(sce)),
  rowData = rd, 
  colData = DataFrame(xyz), 
  spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
  imgData = img,
  sample_id = "sample01"))

# test the plotting
plotVisium(spe, spots=T, point_size=1)

# try the construction from part of the object ----------------------------

# extract the count data
counts <- assay(brain)

# extract the rowData: construct observation & feature metadata
rd2 <- S4Vectors::DataFrame(
  symbol = rowData(brain)$Symbol)

# extract the metadata
xyz2 <- cbind(colData(brain),spatialCoords(brain))

# extract the image
img2 <- imgData(brain)

spe2 <- SpatialExperiment(
  assays = list(counts = counts),
  rowData = rd2, 
  colData = DataFrame(xyz2), 
  spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
  imgData = img2,
  sample_id = "sample01")

# test the plotting
plotVisium(spe2, spots=T, point_size=1)
