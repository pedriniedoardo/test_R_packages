# Introduction to the scDblFinder package ---------------------------------
# Pierre-Luc Germain1* and Aaron Lun**
# 1University and ETH Zürich
# 
# *pierre-luc.germain@hest.ethz.ch
# **infinite.monkeys.with.keyboards@gmail.com
# 
# 15 April 2025
# Abstract
# An introduction to the various methods included in the scDblFinder package.
# 
# Package
# scDblFinder 1.22.0

# 1Introduction -----------------------------------------------------------
# The scDblFinder package gathers various methods for the detection and handling of doublets/multiplets in single-cell sequencing data (i.e. multiple cells captured within the same droplet or reaction volume). This vignette provides a brief overview of the different approaches (which are each covered in their own vignettes) for single-cell RNA sequencing. For doublet detection in genomic data, see the scATACseq vignette. For a more general introduction to the topic of doublets, refer to the OCSA book.
# All methods require as an input either a matrix of counts or a SingleCellExperiment containing count data. With the exception of findDoubletClusters, which operates at the level of clusters (and consequently requires clustering information), all methods try to assign each cell a score indicating its likelihood (broadly understood) of being a doublet.
# The approaches described here are complementary to doublets identified via cell hashes and SNPs in multiplexed samples: while hashing/genotypes can identify doublets formed by cells of the same type (homotypic doublets) from two samples, which are often nearly undistinguishable from real cells transcriptionally (and hence generally unidentifiable through the present package), it cannot identify doublets made by cells of the same sample, even if they are heterotypic (formed by different cell types). Indeed, recent evidence suggests that doublets are for instance a serious and strongly underestimated issue in 10x Flex datasets (see Howitt et al., 2024). Instead, the methods presented here are primarily geared towards the identification of heterotypic doublets, which for most purposes are also the most critical ones.

# 1.1computeDoubletDensity
# The computeDoubletDensity method (formerly scran::doubletCells) generates random artificial doublets from the real cells, and tries to identify cells whose neighborhood has a high local density of articial doublets. See computeDoubletDensity for more information.

# 1.2recoverDoublets
# The recoverDoublets method is meant to be used when some doublets are already known, for instance through genotype-based calls or cell hashing in multiplexed experiments. The function then tries to identify intra-sample doublets that are neighbors to the known inter-sample doublets. See recoverDoublets for more information.

# 1.3scDblFinder
# The scDblFinder method combines both known doublets (if available) and cluster-based artificial doublets to identify doublets. The approach builds and improves on a variety of earlier efforts, and is at present the most accurate approach included in this package. See scDblFinder for more information.

# 1.4directDblClassification
# The directDblClassification method identifies doublets by training a classifier directly on gene expression. This follows the same procedure as scDblFinder for doublet generation and iterative training, but skips the k-nearest neighbor step and directly uses the matrix of real cells and artificial doublets. This is computationally more intensive and generally leads to worse predictions than scDblFinder, and it is included chiefly for comparative purposes. See ?directDblClassification for more information.

# 1.5findDoubletClusters
# The findDoubletClusters method identifies clusters that are likely to be composed of doublets by estimating whether their expression profile lies between two other clusters. See findDoubletClusters for more information.


# 2Installation -----------------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("scDblFinder")
# 
# # or, to get that latest developments:
# BiocManager::install("plger/scDblFinder")

# 3Which method to choose? ------------------------------------------------
# A benchmark of the main methods available in the package is presented in the scDblFinder paper. While the different methods included here have their values, overall the scDblFinder method had the best performance (also superior to other methods not included in this package), and should be used by default.
