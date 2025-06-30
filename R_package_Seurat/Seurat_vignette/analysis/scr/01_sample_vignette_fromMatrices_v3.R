# libraries --------------------------------------------------------------- 
library(Seurat) 
library(cowplot) 
library(harmony) 
library(patchwork) 

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v3")

# Seurat V3 Interface ----------------------------------------------------- 
# This brief vignette demonstrates how to use Harmony with Seurat V3. This example closely follows the Seurat vignette: https://satijalab.org/seurat/v3.0/immune_alignment.html 
# To begin, download the sparse gene count matrices and move them somewhere convenient (e.g. data/). 
load('../data/pbmc_stim.RData') 
ctrl.sparse
stim.sparse

# Initialize Seurat Object ------------------------------------------------ 
# Before running Harmony, make a Seurat object and following the standard pipeline through PCA. 
# IMPORTANT DIFFERENCE: In the Seurat integration tutorial, you need to define a Seurat object for each dataset. With Harmony integration, create only one Seurat object with all cells. 
pbmc <- CreateSeuratObject(counts = cbind(stim.sparse, ctrl.sparse), project = "PBMC", min.cells = 5) %>% 
  Seurat::NormalizeData(verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%  
  ScaleData(verbose = FALSE) %>%  
  RunPCA(pc.genes = pbmc@var.genes, npcs = 20, verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:20) %>%  
  FindNeighbors(reduction = "pca", dims = 1:20) %>%  
  FindClusters(resolution = 0.5) %>%  
  identity() 

# Make sure that the dataset ID is in the object's metadata. Here, we define datasets with the variable stim. 
pbmc@meta.data$stim <- c(rep("STIM", ncol(stim.sparse)), rep("CTRL", ncol(ctrl.sparse))) 
# There is a clear difference between the datasets in the uncorrected PCs 
p1 <- DimPlot(object = pbmc, reduction = "pca", pt.size = .1, group.by = "stim") 
p2 <- VlnPlot(object = pbmc, features = "PC_1", group.by = "stim", pt.size = .1) 
p1 + p2

DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters", pt.size = .1, split.by = 'stim')  

# Run Harmony -------------------------------------------------------------  
# The simplest way to run Harmony is to pass the Seurat object and specify which variable(s) to integrate out. RunHarmony returns a Seurat object, updated with the corrected Harmony coordinates. Let's set plot_convergence to TRUE, so we can make sure that the Harmony objective function gets better with each round.  
pbmc_h <- pbmc %>%   
  RunHarmony("stim", plot_convergence = TRUE)
# Harmony with two or more covariates  
# Do the same with your Seurat object:  
# seuratObject <- RunHarmony(seuratObject, c("dataset", "donor", "batch_id"))  
# To directly access the new Harmony embeddings, use the Embeddings command.  
harmony_embeddings <- Embeddings(pbmc_h, 'harmony')  
harmony_embeddings[1:5, 1:5] 

# Let's make sure that the datasets are well integrated in the first 2 dimensions after Harmony.  
p1 <- DimPlot(object = pbmc_h, reduction = "harmony", pt.size = .1, group.by = "stim")  
p2 <- VlnPlot(object = pbmc_h, features = "harmony_1", group.by = "stim", pt.size = .1)  
p1 + p2  

# Downstream analysis -----------------------------------------------------  
# Many downstream analyses are performed on low dimensional embeddings, not gene expression. To use the corrected Harmony embeddings rather than PCs, set reduction = 'harmony'. For example, let's perform the UMAP and Nearest Neighbor analyses using the Harmony embeddings.  
pbmc_h <- pbmc_h %>%   
  RunUMAP(reduction = "harmony", dims = 1:20) %>%   
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%   
  FindClusters(resolution = 0.5) %>%   
  identity()  
# In the UMAP embedding, we can see more intricate structure. Since we used harmony embeddings, the UMAP embeddings are well mixed.  
DimPlot(pbmc_h, reduction = "umap", group.by = "stim", pt.size = .1, split.by = 'stim')  

DimPlot(pbmc_h, reduction = "umap", label = TRUE, pt.size = .1)  

DimPlot(pbmc_h, reduction = "umap", group.by = "seurat_clusters", pt.size = .1, split.by = 'stim')
