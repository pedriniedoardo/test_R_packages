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
ctrl.sparse[50:70,1:10] 
stim.sparse[50:70,1:10] 

# Initialize Seurat Object ------------------------------------------------ 
ctrl <- CreateSeuratObject(ctrl.sparse,project = "ctrl") %>% 
  Seurat::NormalizeData(verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) 

stim <- CreateSeuratObject(stim.sparse,project = "stim") %>% 
  Seurat::NormalizeData(verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) 

list_obj <- list(ctrl = ctrl, 
                 stim = stim) 

# merge approach ----------------------------------------------------------
# define the ids of the dataset
data.list.id <- names(list_obj)

data.combined.all <- merge(list_obj[[1]], y = list_obj[-1], add.cell.ids = data.list.id, project = "test_harmony")

data.combined.all

# save the sperse matrix
total_sm <- data.combined.all@assays$RNA@counts

# check the total_sm
dim(total_sm)

# confirm the total number of cells
lapply(list_obj,function(x){
  dim(x@assays$RNA@counts)[2]
}) %>%
  unlist() %>%
  sum()

# I need to create a single object to add the cell cycle scoring and other metadata
sobj_total <- CreateSeuratObject(counts = total_sm, project = "PBMC", min.cells = 5)

# add the cell cycle analysis
DefaultAssay(sobj_total) <- "RNA"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sobj_total <- CellCycleScoring(sobj_total, s.features = s.genes, g2m.features = g2m.genes)
sobj_total$percent.mt <- PercentageFeatureSet(sobj_total, pattern = "^MT-")
sobj_total$percent.ribo <- PercentageFeatureSet(sobj_total, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")

# rescale the data for regressing out the sources of variation do not scale all the genes. if needed scale them before the heatmap call
sobj_total <- sobj_total %>%
  Seurat::NormalizeData(verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%  
  ScaleData(verbose = FALSE) %>%  
  RunPCA(pc.genes = pbmc@var.genes, npcs = 20, verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:20) %>%  
  FindNeighbors(reduction = "pca", dims = 1:20) %>%  
  FindClusters(resolution = 0.5) %>%  
  identity()

# Make sure that the dataset ID is in the object's metadata. Here, we define datasets with the variable stim. 
sobj_total@meta.data$stim <-  sobj_total@meta.data$orig.ident

# There is a clear difference between the datasets in the uncorrected PCs 
p1 <- DimPlot(object = sobj_total, reduction = "pca", pt.size = .1, group.by = "stim") 
p2 <- VlnPlot(object = sobj_total, features = "PC_1", group.by = "stim", pt.size = .1) 
p1 + p2 
DimPlot(sobj_total, reduction = "umap", group.by = "seurat_clusters", pt.size = .1, split.by = 'stim') 

# Run Harmony ------------------------------------------------------------- 
# The simplest way to run Harmony is to pass the Seurat object and specify which variable(s) to integrate out. RunHarmony returns a Seurat object, updated with the corrected Harmony coordinates. Let's set plot_convergence to TRUE, so we can make sure that the Harmony objective function gets better with each round. 
pbmc_h <- sobj_total %>%  
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
# In this well mixed embedding, we can start identifying shared cell types with clustering analysis. 
DimPlot(pbmc_h, reduction = "umap", label = TRUE, pt.size = .1) 
DimPlot(pbmc_h, reduction = "umap", group.by = "seurat_clusters", pt.size = .1, split.by = 'stim')
