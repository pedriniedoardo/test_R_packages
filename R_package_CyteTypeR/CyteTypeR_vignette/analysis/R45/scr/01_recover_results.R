# AIM ---------------------------------------------------------------------
# CyteTypeR test recover results

# libraries ---------------------------------------------------------------
library(CyteTypeR)

# read in the object ------------------------------------------------------
# sobj <- readRDS("/idle/ric.cosr/ric.cosr/pedrini.edoardo/sc-rna-seq-cosr-standard-workflow-test2/results/Seurat/object/integrated_obj_cellbender.rds")
token <- Sys.getenv("CYTETYPER_TOKEN")

# processing --------------------------------------------------------------
# this is taken from the log of the failed run 
# https://prod.cytetype.nygen.io/report/b20ef46b-3ba6-428e-b8a6-ac52a3c79f67
# If disconnected, retrieve results with: GetResults()

# pull only the metadata
test_meta <- GetResults(job_id = "b20ef46b-3ba6-428e-b8a6-ac52a3c79f67",
                        auth_token = token)
