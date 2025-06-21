# Introduction ------------------------------------------------------------
# Pau Badia-i-Mompel1 and Jesús Vélez-Santiago2
# 1Heidelberg Universiy
# 2National Autonomous University of Mexico
# 
# 15 April 2025
# Package
# decoupleR 2.14.0
# 
# 1Installation -----------------------------------------------------------
# decoupleR is an R package distributed as part of the Bioconductor project. To install the package, start R and enter:
# install.packages("BiocManager")
# BiocManager::install("decoupleR")
# Alternatively, you can instead install the latest development version from GitHub with:
# BiocManager::install("saezlab/decoupleR")

# 2Usage ------------------------------------------------------------------
# decoupleR (Badia-i-Mompel, Santiago, Braunger, Geiss, Dimitrov, Müller-Dott, Taus, Dugourd, Holland, Flores, and Saez-Rodriguez, 2022) contains different statistical methods to extract biological activities from omics data using prior knowledge. Some of them are:
# AUCell: (Aibar, Bravo Gonzalez-Blas, Moerman, Huynh-Thu, Imrichova, Hulselmans, Rambow, Marine, Geurts, Aerts, van den Oord, Kalender Atak, Wouters, and Aerts, 2017)
# Fast GSEA: (Korotkevich, Sukhov, and Sergushichev, 2019)
# GSVA: (H{ä}nzelmann, Castelo, and Guinney, 2013)
# viper: (Alvarez, Shen, Giorgi, Lachmann, Ding, Ye, and Califano, 2016)
# In this vignette we showcase how to use it with some toy data.

# 2.1Libraries
# decoupleR can be imported as:
library(decoupleR)

# Extra libraries
library(dplyr)
library(pheatmap)
library(fgsea)
library(AUCell)
library(doMC)
library(GSVA)
library(ranger)
library(viper)

# 2.2Input data
# decoupleR needs a matrix (mat) of any molecular readouts (gene expression, logFC, p-values, etc.) and a network that relates target features (genes, proteins, etc.) to “source” biological entities (pathways, transcription factors, molecular processes, etc.). Some methods also require the mode of regulation (MoR) for each interaction, defined as negative or positive weights.
# To get an example data-set, run:
data <- get_toy_data()

mat <- data$mat
head(mat,5)[,1:5]

network <- data$network
network

# This example consists of two small populations of samples (S, cols) with different gene expression patterns (G, rows):
pheatmap(mat, cluster_rows = F, cluster_cols = F)

# Here we can see that some genes seem to be more expressed in one group of samples than in the other and vice-versa. Ideally, we would like to capture these differences in gene programs into interpretable biological entities. In this example we will do it by summarizing gene expression into transcription factor activities.
# The toy data also contains a simple net consisting of 3 transcription factors (Ts) with specific regulation to target genes (either positive or negative). This network can be visualized like a graph. Green edges are positive regulation (activation), red edges are negative regulation (inactivation):
# According to this network, the first population of samples should show high activity for T1 and T3, while the second one only for T2.

# 2.3Methods
# decoupleR contains several methods. To check how many are available, run:
show_methods()

# Each method models biological activities in a different manner, sometimes returning more than one estimate or providing significance of the estimation. To know what each method returns, please check their documentation like this ?run_mlm.

# To have a unified framework, methods have these shared arguments:
# mat : input matrix of molecular readouts.
# network : input prior knowledge information relating molecular features to biological entities.
# .source,.target and .mor : column names where to extract the information from network.
# .source refers to the biological entities.
# .target refers to the molecular features.
# .mor refers to the “strength” of the interaction (if available, else 1s will be used). Only available for methods that can model interaction weights.
# minsize : Minimum of target features per biological entity (5 by default). If less, sources are removed. This filtering prevents obtaining noisy activities from biological entities with very few matching target features in matrix. For this example data-set we will have to keep it to 0 though.

# 2.4Running methods ------------------------------------------------------
# 2.4.1Individual methods
# As an example, let’s first run the Gene Set Enrichment Analysis method (gsea), one of the most well-known statistics:
res_gsea <- run_fgsea(mat, network, .source='source', .target='target', nproc=1, minsize = 0)

res_gsea

# Methods return a result data-frame containing:
# statistic: name of the statistic. Depending on the method, there can be more than one per method.
# source: name of the biological entity.
# condition: sample name.
# score: inferred biological activity.
# p_value: if available, significance of the inferred activity.
# In the case of gsea, it returns a simple estimate of activities (fgsea), a normalized estimate (norm_fgsea) and p-values after doing permutations.
# Other methods can return different things, for example Univariate Linear Model (ulm):
res_ulm <- run_ulm(mat, network, .source='source', .target='target', .mor='mor', minsize = 0)
res_ulm

# In this case, ulm returns just an estimate (ulm) and its associated p-values. Each method can return different statistics, we recommend to check their documentation to know more about them.
# Let us plot the obtained results, first for gsea:

# Transform to matrix
mat_gsea <- res_gsea %>%
  filter(statistic=='fgsea') %>%
  pivot_wider_profile(id_cols = source, names_from = condition, 
                      values_from = score) %>%
  as.matrix()

pheatmap(mat_gsea, cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 40)

# We can observe that for transcription factors T1 and T2, the obtained activities correctly distinguish the two sample populations. T3, on the other hand, should be down for the second population of samples since it is a repressor. This mislabeling of activities happens because gsea cannot model weights when inferring biological activities.
# When weights are available in the prior knowledge, we definitely recommend using any of the methods that take them into account to get better estimates, one example is ulm:
  
# Transform to matrix
mat_ulm <- res_ulm %>%
  filter(statistic=='ulm') %>%
  pivot_wider_profile(id_cols = source, names_from = condition, 
                      values_from = score) %>%
  as.matrix()

pheatmap(mat_ulm, cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 40)

# Since ulm models weights when estimating biological activities, it correctly assigns T3 as inactive in the second population of samples.

# 2.4.2Multiple methods
# decoupleR also allows to run multiple methods at the same time. Moreover, it computes a consensus score based on the obtained activities across methods, called consensus.
# By default, deocuple runs only the top performer methods in our benchmark (mlm, ulm and wsum), and estimates a consensus score across them. Specific arguments to specific methods can be passed using the variable args. For more information check ?decouple.
res_decouple <- decouple(mat, 
                         network, 
                         .source='source', 
                         .target='target',
                         minsize = 0)
res_decouple

# Let us see the result for the consensus score in the previous decouple run:
# Transform to matrix
mat_consensus <- res_decouple %>%
  filter(statistic=='consensus') %>%
  pivot_wider_profile(id_cols = source, names_from = condition, 
                      values_from = score) %>%
  as.matrix()

pheatmap(mat_consensus, cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 40)

# We can observe that the consensus score correctly predicts that T1 and T3 should be active in the first population of samples while T2 in the second one.

# -------------------------------------------------------------------------
# test run all the methods. make sure all the accessory packages are installed
stats <- c("aucell","fgsea","gsva","mdt","mlm","ora","udt","ulm","viper","wmean","wsum")
res_decouple_all <- res_decouple <- decouple(mat,
                                             statistics = stats,
                                             network,
                                             .source='source',
                                             .target='target',
                                             minsize = 0)

# some methods have failed
res_decouple_all %>%
  filter(is.na(score)) %>%
  group_by(statistic) %>%
  summarise()

res_decouple_all %>%
  group_by(statistic) %>%
  summarise()

# Let us see the result for the consensus score in the previous decouple run:
# Transform to matrix
mat_consensus_all <- res_decouple_all %>%
  filter(statistic=='consensus') %>%
  pivot_wider_profile(id_cols = source, names_from = condition, 
                      values_from = score) %>%
  as.matrix()

pheatmap(mat_consensus_all, cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 40)
