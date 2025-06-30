# AIM ---------------------------------------------------------------------
# Integrative analysis in Seurat v5
# Compiled: October 31, 2023
# https://satijalab.org/seurat/articles/seurat5_integration


# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(tidyverse)
library(patchwork)
library(ComplexHeatmap)

options(future.globals.maxSize = 1e9)
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")

# python version ----------------------------------------------------------
Sys.setenv(RETICULATE_PYTHON = "/home/edo/micromamba/envs/env_scvi/bin/python")
library(reticulate)
reticulate::use_condaenv(condaenv = "/home/edo/micromamba/envs/env_scvi/")

py_config()
Sys.getenv()

# -------------------------------------------------------------------------
## Introduction
# Integration of single-cell sequencing datasets, for example across experimental batches, donors, or conditions, is often an important step in scRNA-seq workflows. Integrative analysis can help to match shared cell types and states across datasets, which can boost statistical power, and most importantly, facilitate accurate comparative analysis across datasets. In previous versions of Seurat we introduced methods for integrative analysis, including our ‘anchor-based’ integration workflow. Many labs have also published powerful and pioneering methods, including [Harmony](https://github.com/immunogenomics/harmony) and [scVI](https://yoseflab.github.io/software/scvi-tools/), for integrative analysis. 
# We recognize that while the goal of matching shared cell types across datasets may be important for many problems, users may also be concerned about which method to use, or that integration could result in a loss of biological resolution. In Seurat v5, we introduce more flexible and streamlined infrastructure to run different integration algorithms with a single line of code. This makes it easier to explore the results of different integration methods, and to compare these results to a workflow that excludes integration steps.
# For this vignette, we use a [dataset of human PBMC profiled with seven different technologies](https://www.nature.com/articles/s41587-020-0465-8), profiled as part of a systematic comparative analysis (`pbmcsca`). The data is available as part of our [SeuratData](https://github.com/satijalab/seurat-data) package. 

## Layers in the Seurat v5 object
# Seurat v5 assays store data in layers. These layers can store raw, un-normalized counts (`layer='counts'`), normalized data (`layer='data'`), or z-scored/variance-stabilized data (`layer='scale.data'`). We can load in the data, remove low-quality cells, and obtain predicted cell annotations (which will be useful for assessing integration later), using our [Azimuth pipeline](https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html).

# install the tutorial dataset
# InstallData("pbmcref")
# InstallData("pbmcsca")

# load in the pbmc systematic comparative analysis dataset
obj <- LoadData("pbmcsca")
# notice that no dim reduction is available on the object
# test <- LoadReference("renv/library/linux-ubuntu-jammy/R-4.4/x86_64-pc-linux-gnu/pbmcref.SeuratData/azimuth/")
# DimPlot(test$plot, group.by = "celltype.l2")

# soft filtering
obj  <- subset(obj, nFeature_RNA > 1000)
head(obj@meta.data)

# run the annotation
obj <- RunAzimuth(obj, reference = "pbmcref")
# notice that a dim.reduction is added to the object
DimPlot(obj, group.by = "predicted.celltype.l2")

# currently, the object has two layers in the RNA assay: counts, and data
# notice that after running Azimuth, the object will have a new annotation attached
obj
head(obj@meta.data)

# The object contains data from nine different batches (stored in the `Method` column in the object metadata), representing seven different technologies. We will aim to integrate the different batches together. In previous versions of Seurat, we would require the data to be represented as nine different Seurat objects. When using Seurat v5 assays, we can instead keep all the data in one object, but simply split the layers.
# After splitting, there are now 18 layers (a `counts` and `data` layer for each batch). We can also run a standard scRNA-seq analysis (i.e. without integration). Note that since the data is split into layers, normalization and variable feature identification is performed for each batch independently (a consensus set of variable features is automatically identified).
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Method)
obj

# standard workflow on independent batches
obj <- NormalizeData(obj) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

# We can now visualize the results of a standard analysis without integration. Note that cells are grouping both by cell type and by underlying method. While a UMAP analysis is just a visualization of this, clustering this dataset would return predominantly batch-specific clusters. Especially if previous cell-type annotations were not available, this would make downstream analysis extremely challenging.  
obj <- FindNeighbors(obj, dims=1:30, reduction = 'pca') %>%
  FindClusters(resolution = 2, cluster.name = "unintegrated_clusters") %>%
  RunUMAP(dims = 1:30, reduction = 'pca', reduction.name = 'umap.unintegrated')

# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
DimPlot(obj, reduction = 'umap.unintegrated', group.by=c('Method','predicted.celltype.l2'))

# -------------------------------------------------------------------------
## Perform streamlined (one-line) integrative analysis
# Seurat v5 enables streamlined integrative analysis using the `IntegrateLayers` function. The method currently supports five integration methods. Each of these methods performs integration in low-dimensional space, and returns a dimensional reduction (i.e. `integrated.rpca`) that aims to co-embed shared cell types across batches:

# Anchor-based CCA integration (`method=CCAIntegration`)
# Anchor-based RPCA integration (`method=RPCAIntegration`)
# Harmony (`method=HarmonyIntegration`)
# FastMNN (`method= FastMNNIntegration`)
# scVI (`method=scVIIntegration`)

# Note that our anchor-based RPCA integration represents a faster and more conservative (less correction) method for integration. For interested users, we discuss this method in more detail in our [previous RPCA vignette](https://satijalab.org/seurat/articles/integration_rpca)

# You can find more detail on each method, and any installation prerequisites, in Seurat's documentation (for example, `?scVIIntegration`). For example, scVI integration requires `reticulate` which can be installed from CRAN (`install.packages("reticulate")`) as well as `scvi-tools` and its dependencies installed in a conda environment. Please see scVI installation instructions [here](https://docs.scvi-tools.org/en/stable/installation.html).


# Each of the following lines perform a new integration using a single line of code:
obj <- IntegrateLayers(
  object = obj, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = 'integrated.cca',
  verbose = T)

obj <- IntegrateLayers(
  object = obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = 'integrated.rpca',
  verbose = T)

obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = 'harmony',
  verbose = T)

obj <- IntegrateLayers(
  object = obj, method = FastMNNIntegration,
  new.reduction = 'integrated.mnn',
  verbose = T)

# saveRDS(obj,"../out/integrative_test_obj_V5.rds")

# micromamba create -n env_scvi -c conda-forge r-base r-essentials r-reticulate scvi-tools python scanpy
obj <- IntegrateLayers(
  object = obj, method = scVIIntegration,
  new.reduction = 'integrated.scvi',
  conda_env = '/home/edo/micromamba/envs/env_scvi', verbose = T)

# scvi.reduc <- readRDS("/brahms/haoy/seurat5/object/pbmcsca_scvi.dr.rds")@cell.embeddings
# scvi.reduc <- scvi.reduc[Cells(obj),]
# obj[["integrated.scvi"]] <- CreateDimReducObject(embeddings = scvi.reduc)

# For any of the methods, we can now visualize and cluster the datasets. We show this for CCA integration and Harmony, but you can do this for any method:
obj <- FindNeighbors(obj, reduction = 'integrated.cca', dims = 1:30) %>%
  FindClusters(resolution = 2, cluster.name = 'cca_clusters') %>%
  RunUMAP(reduction = "integrated.cca", dims = 1:30, reduction.name = 'umap.cca')

p1 <- DimPlot(obj, reduction = "umap.cca",
              group.by = c("Method", "predicted.celltype.l2", "cca_clusters"),
              combine = FALSE, label.size = 2) 


obj <- FindNeighbors(obj, reduction = 'harmony', dims = 1:30) %>%
  FindClusters(resolution = 2, cluster.name = 'harmony_clusters') %>%
  RunUMAP(reduction = "harmony", dims = 1:30, reduction.name = 'umap.harmony')

p2 <- DimPlot(
  obj, reduction = "umap.harmony",
  group.by = c("Method", "predicted.celltype.l2", "harmony_clusters"),
  combine = FALSE, label.size = 2)

wrap_plots(c(p1, p2), ncol = 2, byrow = F)

obj <- FindNeighbors(obj, reduction = 'integrated.rpca', dims = 1:30) %>%
  FindClusters(resolution = 2, cluster.name = 'rpca_clusters') %>%
  RunUMAP(reduction = "integrated.rpca", dims = 1:30, reduction.name = 'umap.rpca')


# We hope that by simplifying the process of performing integrative analysis, users can more carefully evaluate the biological information retained in the integrated dataset. For example, users can compare the expression of biological markers based on different clustering solutions, or visualize one method's clustering solution on different UMAP visualizations.

p12 <- VlnPlot(obj, features = "rna_CD8A", group.by = 'unintegrated_clusters') + NoLegend() + ggtitle("CD8A - Unintegrated Clusters")
p22 <- VlnPlot(obj, "rna_CD8A", group.by = 'cca_clusters') + NoLegend() + ggtitle("CD8A - CCA Clusters")
p32 <- VlnPlot(obj, "rna_CD8A", group.by = 'harmony_clusters') + NoLegend() + ggtitle("CD8A - Harmony Clusters")

p12 | p22 | p32

p4 <- DimPlot(obj, reduction="umap.unintegrated", group.by=c("cca_clusters"))
p5 <- DimPlot(obj, reduction="umap.rpca", group.by=c("cca_clusters"))
p6 <- DimPlot(obj, reduction="umap.harmony", group.by=c("cca_clusters"))
p4 | p5 | p6

# see the clustering with the cell annotation from azimuth
p42 <- DimPlot(obj, reduction="umap.unintegrated", group.by=c("predicted.celltype.l1"))
p52 <- DimPlot(obj, reduction="umap.rpca", group.by=c("predicted.celltype.l1"))
p62 <- DimPlot(obj, reduction="umap.harmony", group.by=c("predicted.celltype.l1"))
p42 | p52 | p62


# Once integrative analysis is complete, you can rejoin the layers - which collapses the individual datasets together and recreates the original `counts` and `data` layers. You will need to do this before performing any differential expression analysis. However, you can always resplit the layers in case you would like to reperform integrative analysis.
obj <- JoinLayers(obj)
obj

# -------------------------------------------------------------------------
# check cross annotation
df_harmony.l1 <- obj@meta.data %>%
  group_by(harmony_clusters,predicted.celltype.l1) %>%
  summarise(n = n(), .groups = 'drop') %>%
  # add entries with 0 occurrencies
  complete(harmony_clusters,predicted.celltype.l1, fill = list(n = 0)) %>%
  group_by(harmony_clusters) %>%
  mutate(tot_cluster = sum(n)) %>%
  ungroup() %>%
  mutate(prop = n/tot_cluster)

# plot the data
mat_harmony.l1 <- df_harmony.l1 %>%
  select(harmony_clusters,predicted.celltype.l1,prop) %>%
  pivot_wider(names_from = predicted.celltype.l1, values_from = prop) %>%
  column_to_rownames("harmony_clusters")

Heatmap(mat_harmony.l1)

# suggests the most common annotation per cluster
df_harmony.l1_top <- df_harmony.l1 %>%
  group_by(harmony_clusters) %>%
  slice_max(prop,n = 1,with_ties = F) %>%
  ungroup()

df_meta <- obj@meta.data

obj$harmony_annotation <- df_meta %>%
  left_join(df_harmony.l1_top %>% select(harmony_clusters,harmony_annotation = predicted.celltype.l1), by = "harmony_clusters") %>%
  pull(harmony_annotation)

p61 <- DimPlot(obj, reduction="umap.harmony", group.by=c("harmony_clusters"))
p62 <- DimPlot(obj, reduction="umap.harmony", group.by=c("predicted.celltype.l1"))
p63 <- DimPlot(obj, reduction="umap.harmony", group.by=c("harmony_annotation"))
p61 | p62 | p63

df_cca.l1 <- obj@meta.data %>%
  group_by(cca_clusters,predicted.celltype.l1) %>%
  summarise(n = n(), .groups = 'drop') %>%
  # add entries with 0 occurrencies
  complete(cca_clusters,predicted.celltype.l1, fill = list(n = 0)) %>%
  group_by(cca_clusters) %>%
  mutate(tot_cluster = sum(n)) %>%
  ungroup() %>%
  mutate(prop = n/tot_cluster)

# plot the data
mat_cca.l1 <- df_cca.l1 %>%
  select(cca_clusters,predicted.celltype.l1,prop) %>%
  pivot_wider(names_from = predicted.celltype.l1, values_from = prop) %>%
  column_to_rownames("cca_clusters")

Heatmap(mat_cca.l1)

# suggests the most common annotation per cluster
df_cca.l1_top <- df_cca.l1 %>%
  group_by(cca_clusters) %>%
  slice_max(prop,n = 1,with_ties = F) %>%
  ungroup()

obj$cca_annotation <- df_meta %>%
  left_join(df_cca.l1_top %>% select(cca_clusters,cca_annotation = predicted.celltype.l1), by = "cca_clusters") %>%
  pull(cca_annotation)

p71 <- DimPlot(obj, reduction="umap.cca", group.by=c("cca_clusters"))
p72 <- DimPlot(obj, reduction="umap.cca", group.by=c("predicted.celltype.l1"))
p73 <- DimPlot(obj, reduction="umap.cca", group.by=c("cca_annotation"))
p71 | p72 | p73

# -------------------------------------------------------------------------
# # Identify conserved cell type markers
# # To identify canonical cell type marker genes that are conserved across conditions, we provide the FindConservedMarkers() function. This function performs differential gene expression testing for each dataset/group and combines the p-values using meta-analysis methods from the MetaDE R package. For example, we can calculated the genes that are conserved markers irrespective of stimulation condition in cluster 6 (NK cells).
# Idents(obj) <- "harmony_annotation"
# CD8.markers <- FindConservedMarkers(obj, ident.1 = "CD8 T", grouping.var = "Method", verbose = T)
# head(CD8.markers)
# # You can perform these same analysis on the unsupervised clustering results (stored in seurat_clusters), and use these conserved markers to annotate cell types in your dataset.
# # The DotPlot() function with the split.by parameter can be useful for viewing conserved cell type markers across conditions, showing both the expression level and the percentage of cells in a cluster expressing any given gene. Here we plot 2-3 strong marker genes for each of our 14 clusters.
# 
# color_id <- alphabet(length(unique(obj@meta.data$Method)))
# names(color_id) <- unique(obj@meta.data$Method)
# 
# # check the colors
# show_col(color_id)
# 
# DotPlot(obj,idents = c("NK","CD8 T"),
#         features = c("CCL5","CD8A","GZMH","NKG7","CST7","CD8B","B2M","HLA-C","FGFBP2","CD3D"),
#         cols = color_id,
#         dot.scale = 8,
#         split.by = "Method") +
#   RotatedAxis()

Idents(obj) <- "harmony_annotation"
harmony_markers_genes <- FindAllMarkers(obj, verbose = T)

Idents(obj) <- "cca_annotation"
cca_markers_genes <- FindAllMarkers(obj, verbose = T)

markers_gene <- full_join(harmony_markers_genes %>% select(gene,cluster,avg_log2FC,p_val_adj),
                          cca_markers_genes %>% select(gene,cluster,avg_log2FC,p_val_adj),
                          by = c("gene","cluster"),
                          suffix = c("_harmony","_cca")) %>%
  arrange(cluster)

markers_gene %>%
  ggplot(aes(x=avg_log2FC_harmony,y=avg_log2FC_cca)) +
  geom_point(shape = 1,alpha = 0.5) +
  facet_wrap(~cluster) +
  theme_bw() +
  theme(strip.background = element_blank())+
  geom_abline(intercept = 0,slope = 1,linetype = 2,col="red")

# -------------------------------------------------------------------------
# Lastly, users can also perform integration using sctransform-normalized data (see our [SCTransform vignette](https://satijalab.org/seurat/articles/sctransform_vignette) for more information), by first running SCTransform normalization, and then setting the `normalization.method` argument in `IntegrateLayers`.
obj2 <- LoadData("pbmcsca")
obj2 <- subset(obj2, nFeature_RNA > 1000)
obj2[["RNA"]] <- split(obj2[["RNA"]], f = obj2$Method)

options(future.globals.maxSize = 3e+09)
obj2 <- SCTransform(obj2)
obj2 <- RunPCA(obj2, npcs = 30, verbose = F)
obj2 <- IntegrateLayers(object = obj2, 
                        method = RPCAIntegration,
                        normalization.method="SCT",
                        verbose = F)
obj2 <- FindNeighbors(obj2, dims = 1:30,reduction = 'integrated.dr') %>%
  FindClusters(resolution = 2) %>%
  RunUMAP(dims = 1:30,reduction = 'integrated.dr', reduction.name = 'umap.dr')

DimPlot(obj2, reduction="umap.dr", group.by=c("SCT_snn_res.2"))


# -------------------------------------------------------------------------
# perform differential expression
obj2 <- PrepSCTFindMarkers(obj2)
# obj2$celltype.stim <- paste(obj2$seurat_annotations, obj2$stim, sep = "_")
Idents(obj2) <- "SCT_snn_res.2"
markers_genes <- FindAllMarkers(obj2, verbose = T)

