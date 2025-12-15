# AIM ---------------------------------------------------------------------
# test the tool using a visium spatial dataset 
# this is a sample routine for a seurat object
# here I am loading the file from a regular output of spaceranger

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)

# try to load a full object -----------------------------------------------
# identify the location of the output matrix
matrix_dir <- "../data/sample_Visium/results/spaceranger/merged/V1_Mouse_Brain_Sagittal_Anterior_Section_1/outs/filtered_feature_bc_matrix/"

# read in the matrix data
counts <- Seurat::Read10X(data.dir = matrix_dir)  

# create the object
data <- Seurat::CreateSeuratObject(
  counts = counts , 
  project = 'test', 
  assay = 'Spatial')

# add metadata
data$slice <- 1 
data$region <- 'test' 

# locate image information
imgpath <- "../data/sample_Visium/results/spaceranger/merged/V1_Mouse_Brain_Sagittal_Anterior_Section_1/outs/spatial"
# read in the image information
img <- Seurat::Read10X_Image(image.dir = imgpath)  

# add the image slot
Seurat::DefaultAssay(object = img) <- 'Spatial'  

# add the coldata information
img <- img[colnames(x = data)]
data[['image']] <- img  

# sample plotting
SpatialFeaturePlot(data, features = "nCount_Spatial")
