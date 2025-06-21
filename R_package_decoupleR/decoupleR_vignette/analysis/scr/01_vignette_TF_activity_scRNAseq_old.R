# https://www.bioconductor.org/packages/release/bioc/vignettes/decoupleR/inst/doc/tf_sc.html
# Transcription factor activity inference from scRNA-seq
# Pau Badia-i-Mompel1

# scRNA-seq yield many molecular readouts that are hard to interpret by themselves. One way of summarizing this information is by inferring transcription factor (TF) activities from prior knowledge.

# In this notebook we showcase how to use decoupleR for transcription factor activity inference with a down-sampled PBMCs 10X data-set. The data consists of 160 PBMCs from a Healthy Donor. The original data is freely available from 10x Genomics here from this webpage.
# https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

# 1Loading packages -------------------------------------------------------
# First, we need to load the relevant packages, Seurat to handle scRNA-seq data and decoupleR to use statistical methods.

## We load the required packages
library(Seurat)
library(decoupleR)
library(tidyverse)

# # Only needed for data handling and plotting
# library(dplyr)
# library(tibble)
# library(tidyr)
# library(patchwork)
# library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)

# 2Loading the data-set ---------------------------------------------------
# Here we used a down-sampled version of the data used in the Seurat vignette. We can open the data like this:
# inputs_dir <- system.file("extdata", package = "decoupleR")
# data <- readRDS(file.path(inputs_dir, "sc_data.rds"))
# save it in local 
# saveRDS(data,file = "data/data.rds")
# this is a preprocessed seurat object
data <- readRDS(file = "../data/data.rds")

# We can observe that we have different cell types:
DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# 3DoRothEA network -------------------------------------------------------
# DoRothEA is a comprehensive resource containing a curated collection of TFs and their transcriptional targets. Since these regulons were gathered from different types of evidence, interactions in DoRothEA are classified in different confidence levels, ranging from A (highest confidence) to D (lowest confidence). 
# Moreover, each interaction is weighted by its confidence level and the sign of its mode of regulation (activation or inhibition). The direction is coded in the mor column of the dataset.

# For this example we will use the human version (mouse is also available) and we will use the confidence levels ABC. We can use decoupleR to retrieve it from OmniPath:
net <- get_dorothea(organism='human', levels=c('A', 'B', 'C'))
net

# 4Activity inference with Weighted Mean ----------------------------------
# To infer activities we will run the Weighted Mean method (wmean). It infers regulator activities by first multiplying each target feature by its associated weight which then are summed to an enrichment score wmean. Furthermore, permutations of random target features can be performed to obtain a null distribution that can be used to compute a z-score norm_wmean, or a corrected estimate corr_wmean by multiplying wmean by the minus log10 of the obtained empirical p-value.

# In this example we use wmean but we could have used any other. To see what methods are available use show_methods().

# To run decoupleR methods, we need an input matrix (mat), an input prior knowledge network/resource (net), and the name of the columns of net that we want to use.

# Extract the normalized log-transformed counts
mat <- as.matrix(data@assays$RNA@data)

# Run wmean
acts <- run_wmean(mat=mat, net=net, .source='source', .target='target',
                  .mor='mor', times = 100, minsize = 5)
acts

# 5Visualization ----------------------------------------------------------
# From the obtained results, we will select the norm_wmean activities and store them in our object as a new assay called tfswmean:

# Extract norm_wmean and store it in tfswmean in pbmc
data[['tfswmean']] <- acts %>%
  filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = data) <- "tfswmean"

# Scale the data
data <- ScaleData(data)
data@assays$tfswmean@data <- data@assays$tfswmean@scale.data

# This new assay can be used to plot activities. Here we observe the activity inferred for PAX5 across cells, which it is particulary active in B cells. Interestingly, PAX5 is a known TF crucial for B cell identity and function. The inference of activities from “foot-prints” of target genes is more informative than just looking at the molecular readouts of a given TF, as an example here is the gene expression of PAX5, which is not very informative by itself:
DefaultAssay(object = data)
p1 <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend() + ggtitle('Cell types')
p2 <- (FeaturePlot(data, features = c("PAX5")) & scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('PAX5 activity')

DefaultAssay(object = data) <- "RNA"
p3 <- FeaturePlot(data, features = c("PAX5")) + ggtitle('PAX5 expression')
DefaultAssay(object = data) <- "tfswmean"
p1 | p2 | p3

# 6Exploration ------------------------------------------------------------
# We can also see what is the mean activity per group of the top 20 more variable TFs:
n_tfs <- 25

# Extract activities from object as a long dataframe
df <- t(as.matrix(data@assays$tfswmean@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(data)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
tfs <- df %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
palette_length = 100
my_color <- colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot
# pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 
# implement it in complexheatmap
Heatmap(top_acts_mat)

# Here we can observe other known marker TFs appearing, PAX5 for B cells RFX5 for the myeloid lineage and JUND for the lymphoid.