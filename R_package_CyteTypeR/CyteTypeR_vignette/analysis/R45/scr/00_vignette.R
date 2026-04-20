# AIM ---------------------------------------------------------------------
# run the sample vigentte

# libraries ---------------------------------------------------------------
library(CyteTypeR)
library(dplyr)
library(patchwork)
library(Matrix)
library(Seurat)
library(presto)
library(duckdb)

# Quick Start -------------------------------------------------------------

# prepped_data <- PrepareCyteTypeR(pbmc,
#                                  pbmc.markers,
#                                  n_top_genes = 10,
#                                  group_key = 'seurat_clusters',
#                                  aggregate_metadata = TRUE,
#                                  coordinates_key = "umap")
# 
# metadata <- list(
#   title = 'My scRNA-seq analysis of human pbmc',
#   run_label = 'initial_analysis',
#   experiment_name = 'pbmc_human_samples_study')
# 
# results <- CyteTypeR(prepped_data = prepped_data, 
#                      study_context = "pbmc blood samples from humans", 
#                      metadata = metadata
# )

# configuration -----------------------------------------------------------
# Use a named list for llm configs
# llm_configs=list(
#   "provider" = "openai",   # one of: anthropic, bedrock, google, groq, mistral, openai, openrouter
#   "name" = "gpt-4o-mini",
#   "apiKey" = "your-api-key",
#   "baseUrl" = "https://api.openai.com/v1",   # optional
#   "modelSettings" = list(                       # optional
#     "temperature" = 0.0,
#     "max_tokens" = 4096
#   )
# )
# 
# 
# result <- CyteTypeR(obj = pbmc, 
#                     prepped_data = prepped_data,
#                     study_context = "pbmc blood samples from humans",
#                     metadata = metadata,
#                     llm_configs = llm_configs
# )

# read in a sample dataset ------------------------------------------------
# Load the dataset
pbmc.data <- Read10X(data.dir = "../../data/pbmc3k/filtered_gene_bc_matrices/hg19/")

# Pre-processing ----------------------------------------------------------
# Current version of CyteTypeR works with Seurat objects and requires minimally some basic pre-processing before CyteTypeR can be used.

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200) %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

# all.genes <- rownames(pbmc)
# pbmc <- ScaleData(pbmc, features = all.genes)
# do not scale all the data
pbmc <- ScaleData(pbmc)

# Cluster the cells and run UMAP
pbmc <- RunPCA(pbmc) %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:10)

# Find markers for all Clusters
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)

# save the table for the analysis
pbmc.markers %>%
  write_tsv(file = "../../out/R45/table/pbmc_markers.tsv")

# check the top poritive markers per cluster
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

pbmc@meta.data

DimPlot(pbmc,label = T)

# save the output object processed
saveRDS(pbmc,"../../out/R45/object/pbmc.rds")

# Run CyteTypeR
# Prep data for job submission to cytetype api

# obj: A Seurat object with normalized data and clustering results.
# marker_table: Data frame containing marker genes with columns: 'cluster', 'gene', and 'avg_log2FC'. Output from Seurat's FindAllMarkers() function.
# group_key: Character string specifying the metadata column containing cluster assignments. Default is "seurat_clusters".
# gene_symbols: Character string specifying the gene symbol field name. Default is 'gene_symbols'.
# n_top_genes: Integer specifying the maximum number of top marker genes per cluster to include (filtered by avg_log2FC > 1). Default is 50.
# aggregate_metadata" Logical indicating whether to aggregate metadata across cells within each cluster. Default is TRUE.
# coordinates_key: Character string specifying which dimensional reduction to use for visualization coordinates (e.g., "umap", "tsne"). Default is "umap".
# vars_h5_path: Character string specifying the local file path for the generated vars.h5 artifact (feature expression). Default is "vars.h5".
# obs_duckdb_path: Character string specifying the local file path for the generated obs.duckdb artifact (cell metadata). Default is "obs.duckdb".
prepped_data <- PrepareCyteTypeR(obj = pbmc,
                                 marker_table = pbmc.markers,
                                 n_top_genes = 10,
                                 group_key = 'seurat_clusters',
                                 aggregate_metadata = TRUE,
                                 coordinates_key = "umap",
                                 vars_h5_path = "../../out/R45/object/vars.h5",
                                 obs_duckdb_path = "../../out/R45/object/obs.duckdb")

# Adding metadata on the report
metadata <- list(
  title = 'My scRNA-seq analysis of human pbmc',
  run_label = 'initial_analysis',
  experiment_name = 'pbmc_human_samples_study')


# Submit job to cytetype
# obj: A Seurat object (will be returned with annotations added).
# prepped_data: Named list containing prepared data from PrepareCyteTypeR().
# study_context: Optional character. Biological context for the experimental setup (e.g. organisms, tissues, diseases, developmental stages, single_cell methods, experimental conditions). Default is NULL.
# llm_configs: Optional list of LLM configs. Each element must match the LLMModelConfig schema with required provider and name; either apiKey or all AWS credentials (awsAccessKeyId, awsSecretAccessKey, awsDefaultRegion) must be provided. Default is NULL (API default model).
# n_parallel_clusters: Integer. Number of parallel requests to the model (max 50). High values can trigger rate limits. Default is 2.
# save_query: Logical. Whether to save the request payload to a JSON file. Default is TRUE.
# query_filename: Character. Filename for the saved query when save_query is TRUE. Default is "query.json".
# override_existing_results: Logical. If TRUE, allow overwriting existing results with the same results_prefix. If FALSE and results exist, an error is raised. Default is FALSE.
# require_artifacts: Logical. If TRUE, an error during artifact build or upload stops the run; if FALSE, failures are skipped and annotation continues without artifacts. Default is TRUE.
results <- CyteTypeR(obj = pbmc,
                     prepped_data = prepped_data, 
                     study_context = "pbmc blood samples from humans", 
                     metadata = metadata,
                     query_filename = "../../out/R45/object/query.json")

results@misc
