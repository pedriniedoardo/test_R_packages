# AIM ---------------------------------------------------------------------
# test the tool using a visium spatial dataset

# libraries ---------------------------------------------------------------
library(Seurat)
library(BayesPrism)
library(tidyverse)
library(SeuratData)
library(AnnotationHub)

# read in the data --------------------------------------------------------
# read in the reference sc dataset
seurat_ref <- readRDS("../data/allen_cortex.rds")

# explore the metadata
seurat_ref@meta.data

# plot the dataset with the relative annotaiton
# DimPlot(seurat_ref,group.by = "class")

# read in the sample spatial dataset
# options(timeout=600)
# InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")
# plot to check the object
SpatialFeaturePlot(brain, features = "nCount_Spatial")

# to speed up the computation select a small subset of the object
# add the coordinates
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


# use each barcode as a individual bulk samples for the deconvolution
bulk_counts_full <- GetAssayData(brain,slot = "count") %>%
  # coerce to matrix
  as.matrix()

# to speed-up the computation reduce the number of spots
bulk_counts <- GetAssayData(brain_small,slot = "count") %>%
  # coerce to matrix
  as.matrix()

# this ia a dataset of 2696 samples
dim(bulk_counts)

LUT_prop <- seurat_ref@meta.data %>%
  group_by(class) %>%
  summarise(n = n(),.groups = "drop") %>%
  # group_by(stim) %>%
  mutate(tot = sum(n),
         prop = n/tot)

# Check the dimensions of our simulated data
dim(seurat_ref)
dim(bulk_counts)

#  Step 2: Prepare Data for BayesPrism ------------------------------------
# BayesPrism requires the data in specific formats. We need to extract the raw counts and cell labels from our Seurat object and ensure the gene names match between our reference and bulk data.

# Extract Information from the Seurat Object
# BayesPrism works best with raw counts.
sc_counts_matrix <- GetAssayData(seurat_ref, slot = "counts")

# Get cell type labels. These will be our 'cell.type.labels'.
cell_types <- seurat_ref$class

# Get cell state labels. For a simple case, these can be the same as cell types.
# In more complex analyses, these could be cluster IDs or subtypes.
# https://github.com/Danko-Lab/BayesPrism/issues/66
cell_states <- seurat_ref$class

# Filter and Align Genes
# reshape the data to be accepted by the tool
# generally data are reporeted as feature X sample, but in the tool they are handled as sample X feature.
# we can input sparse matrices: https://github.com/Danko-Lab/BayesPrism/issues/58
sc.dat <- t(sc_counts_matrix)
sc.dat[1:5,1:5]
dim(sc.dat)

bk.dat <- t(bulk_counts)
bk.dat[,1:5]
dim(bk.dat)

# Run the reccommended filtering on the genes in this step for a real data analysis
# Filter outlier genes from scRNA-seq data
# Next, we remove the genes from selected groups. Note that when sex is not identical between the reference and mixture, we recommend excluding genes from chrX and chrY. We also remove lowly transcribed genes, as the measurement of transcription of these genes tend to be noise-prone. Removal of these genes can also speed up computation.
# if using mm do not include MALAT1 category: https://github.com/Danko-Lab/BayesPrism/issues/89
sc.dat.filtered <- cleanup.genes(input = sc.dat,
                                 input.type = "count.matrix",
                                 species = "mm",
                                 # gene.group = c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                 gene.group = c( "Rb","Mrp","other_Rb","chrM","chrX","chrY") ,
                                 exp.cells = 5)

dim(sc.dat.filtered)
dim(sc.dat)

# Next, we check the concordance of gene expression for different types of genes. As bulk and single cell data are usually collected by different experimental protocols, they may have different sensitivity to different types of genes.
# note this function only works for human data. For other species, you are advised to make plots by yourself.

# check if there are genes in common
intersect(colnames(sc.dat.filtered),colnames(bk.dat))

# Note that this function only works for human data (using GENCODE v22 to match TCGA). For other species, you are advised to make plots by yourself.
# plot.bulk.vs.sc(sc.input = sc.dat.filtered,
#                 bulk.input = bk.dat
#                 #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
# )

# Subset protein coding genes.
# A function to select genes by gene type, as informed by the plot from plot.bulk.vs.sc. Only works for human data, as the annotation was for TCGA (GENCODE v22). For other species please filter manually if needed.
# sc.dat.filtered.pc <-  select.gene.type(sc.dat.filtered,
#                                         gene.type = "protein_coding")
# pull the mouse gene annotation
ah <- AnnotationHub()

# unique(ah$species) %>% unique() %>% str_subset(pattern = "uscul")
mouse_ens <- query(ah, c("Mus musculus", "EnsDb"))
mouse_ens <- mouse_ens[["AH119358"]]

LUT_genes <- genes(mouse_ens, return.type = "data.frame")  %>%
  dplyr::select(gene_id, gene_name, entrezid, gene_biotype)

id_gene_filter <- LUT_genes %>%
  dplyr::filter(gene_biotype %in% c("protein_coding")) %>%
  pull(gene_name)

# filter the genes
sc.dat.filtered.pc <- sc.dat.filtered[,intersect(id_gene_filter,colnames(sc.dat.filtered))]

# compare the dimensions
str(sc.dat.filtered.pc)
dim(sc.dat.filtered.pc)
dim(sc.dat.filtered)
dim(sc.dat)

# potentially run the DE to select markers gene per cell_type

# Step 3: Run BayesPrism Deconvolution ------------------------------------
# With the data properly formatted and aligned, we can now run the deconvolution.

# Construct the BayesPrism Object
# We combine the reference and cell labels into a 'prism' object.
my_prism <- new.prism(
  reference = sc.dat.filtered.pc,
  mixture = bk.dat,
  key = NULL, # No specific tumor/normal key needed for this example
  cell.type.labels = cell_types,
  cell.state.labels = cell_states,
  input.type = "count.matrix"
)

# Run the Deconvolution Algorithm
# This is the core computational step.
# We'll use 4 cores for this example. Adjust 'n.cores' as needed.
bp_results <- run.prism(prism = my_prism, n.cores = 20)

# You can look at the returned object to see what it contains
bp_results

# saveRDS(bp_results,"../out/object/bp_results_rawCount_pBulk_treat_100_ifnb.rds")
saveRDS(bp_results,"../out/object/bp_results_test_visium.rds")

# Step 4: Extract and Calculate Deconvoluted Reads ------------------------
# this can be extracted from the object as follows
# x <- "Endothelial"

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

# notice that the number of genes matches the one present in the reference
lapply(list_bulk, function(x){
  dim(x)
})

dim(sc.dat.filtered.pc)

# get the fraction estimated
theta <- get.fraction(bp=bp_results,
                      which.theta="final",
                      state.or.type="type")
theta

# test --------------------------------------------------------------------
# load the matrices in a seurat object
SpatialFeaturePlot(brain_small, features = "nCount_Spatial") + theme(legend.position = "right")

lapply(list_bulk_integer,function(x){
  mat <- x %>%
    as.sparse()
})

# define the paths in input and output
path_in_image <- paste0("../data/misc/test_edo/raw/sample_",samp,"/",samp,"_spatial/")
path_in_reads <- paste0("../data/misc/test_edo/raw/sample_",samp,"/")
sample_h5 <- paste0(samp,".h5")
path_in_meta <- paste0("../data/misc/test_edo/raw/sample_",samp,"/",samp,"_metadata.csv")

# read in the image file. load the highres image
img <- Read10X_Image(image.dir = path_in_image,
                     image.name = "tissue_hires_image.png")

# load the full matrix of reads
spobj <- Load10X_Spatial(path_in_reads,
                         filename = sample_h5,
                         image = img,
                         filter.matrix = F)