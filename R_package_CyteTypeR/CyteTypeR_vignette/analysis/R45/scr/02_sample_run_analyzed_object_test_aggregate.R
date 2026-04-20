# AIM ---------------------------------------------------------------------
# the try to run the tool using a pre-analyzed object
# test the effect of running aggregate_metadata = FALSE

# libraries ---------------------------------------------------------------
library(CyteTypeR)
library(tidyverse)
library(patchwork)
library(Matrix)
library(Seurat)
library(presto)
library(duckdb)

# parameters --------------------------------------------------------------
# input:
rds <- "../../out/R45/object/pbmc.rds"
markers <- "../../out/R45/table/pbmc_markers.tsv"
metadata <- "../../data/meta_annotation_pbmc.json"

message("input rds: ", rds)
message("input markers: ", rds)
# message("input metadata: ", metadata)

# output:
out_id_rds <- "../../out/R45/object/pbmc_CyteTypeR_testAggrFalse.rds"
out_id_varsh5 <- "../../out/R45/object/vars_testAggrFalse.h5"
out_id_duckdb <- "../../out/R45/object/obs_testAggrFalse.duckdb"
out_id_filename <- "../../out/R45/object/query_testAggrFalse.json"

message("output CyteTypeR rds: ", out_id_rds)
message("output CyteTypeR vars.h5: ", out_id_varsh5)
message("output CyteTypeR obs.duckdb: ", out_id_duckdb)
message("output CyteTypeR query.json: ", out_id_filename)

# params
# token <- snakemake@params$token
token <- Sys.getenv("CYTETYPER_TOKEN")
# this should be handled in the snakemake rule
if (token == "") {
  stop("CYTETYPER_TOKEN environment variable is not set or is empty. Please export it before running Snakemake.")
}

# define the test covariate for the sample run
cov_markers <- "seurat_clusters"

# read input --------------------------------------------------------------
# read in the full object
scobj <- readRDS(rds)

# read in the makrers
df_markers <- read_tsv(markers)

# read in the metadata
list_metadata <- yaml::read_yaml(metadata)

# ======================================================================
# == run the standard processing ==
# ======================================================================

# reshape the marker output to be accepted
# This one assumes that the order of the genes is maintained from the output of Seurat's FindAllMarkers() function during the integration step
df_markers_fix <- df_markers %>%
  # mutate(cluster = as.numeric(cluster)) %>%
  # filter(resolution == str_remove(cov_markers,pattern = "RNA_snn_res.")) %>%
  select(cluster,gene,avg_log2FC)

# run the query preparation
prepped_data <- PrepareCyteTypeR(obj = scobj,
                                 marker_table = df_markers_fix,
                                 n_top_genes = 10,
                                 group_key = cov_markers,
                                 aggregate_metadata = FALSE,
                                 coordinates_key = "umap",
                                 vars_h5_path = out_id_varsh5,
                                 obs_duckdb_path = out_id_duckdb)

# run the annotation
results <- CyteTypeR(obj = scobj,
                     auth_token = token,
                     prepped_data = prepped_data, 
                     study_context = list_metadata$study_context, 
                     metadata = list_metadata,
                     query_filename = out_id_filename)

# ======================================================================
# == save output ==
# ======================================================================

saveRDS(results,out_id_rds)
