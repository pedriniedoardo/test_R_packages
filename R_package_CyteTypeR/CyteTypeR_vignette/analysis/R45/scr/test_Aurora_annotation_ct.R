library("Seurat")
library("openxlsx")
library("CyteTypeR")
library("dplyr")
setwd("/idle/ric.cosr/ric.cosr/maurizio.aurora/bop/scrnaseq_cosr_standard_workflow_cr10/results/Seurat/object")
getwd()


integrated_obj <-readRDS("integrated_obj.rds")
DimPlot(integrated_obj, group.by = "orig.ident")
DimPlot(integrated_obj, group.by = "RNA_snn_res.0.3")
Idents(integrated_obj) <- integrated_obj$RNA_snn_res.0.3
DefaultAssay(integrated_obj) <- "RNA"
integrated_obj.markers <- FindAllMarkers(integrated_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


filename_xls <- 'FindAllMarkers_RNA_snn_res.0.3_cb.xlsx'
write.xlsx(integrated_obj.markers,
           file= filename_xls,
           rowNames = T,
           asTable = T)

prepped_data <- PrepareCyteTypeR(
  integrated_obj,
  integrated_obj.markers,
  n_top_genes = 30,
  group_key = 'RNA_snn_res.0.3',
  aggregate_metadata = TRUE,
  coordinates_key = "umap"
)

# questa parte si può migliorare, aggiungendo nel metadata dell' oggetto
# in input una colonna Condition/Contrast etc.
# e poi specificarla qua tipo condition_key = 'Condition'
# devo ancora provare

metadata <- list(
  title = 'scRNA-seq of human LN slices',
  run_label = 'RNA_snn_res.0.3',
  experiment_name = 'Fergusson et al. reanalysis'
)

results <- CyteTypeR(
  obj=integrated_obj,
  prepped_data = prepped_data,
  require_artifacts = FALSE,
  study_context = "human LN slices",
  metadata = metadata
)


saveRDS(results, "CyteTypeR_annotation.Rds")
# seems the best
DimPlot(results, group.by = "cytetype_cellOntologyTerm_RNA_snn_res.0.3")

results$cytetype_annotation_RNA_snn_res.0.3



