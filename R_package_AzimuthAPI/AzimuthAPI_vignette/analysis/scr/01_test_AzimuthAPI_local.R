# AzimuthAPI - local test (ANNotate via reticulate + panhumanpy) ------------
# Mirrors the "Pan-human Azimuth in R" vignette, Option 2 (local run),
# on the bmcite demo dataset (Stuart, Butler et al Cell 2019) from SeuratData.
# https://github.com/satijalab/AzimuthAPI
#
# Assumes packages from 00_libraries.R are already installed via renv, and
# that the `pan-human-azimuth` conda/python environment (see 00_libraries.R)
# already exists with panhumanpy installed.
#
# Run with the working directory set to this project's `analysis/` folder
# (e.g. open analysis.Rproj in RStudio).

library(AzimuthAPI)
library(Seurat)
library(SeuratData)
library(reticulate)
library(tidyverse)
library(presto)


# parameters --------------------------------------------------------------
out_obj_dir <- "../out/object/"
out_plot_dir <- "../out/plot/"

# 1. Load demo data ----------------------------------------------------------
# One-time download, only needed the first time (skip if already installed):
# InstallData("bmcite")

bmcite <- LoadData("bmcite")
bmcite <- NormalizeData(bmcite)

# 2. Point reticulate at the local panhumanpy Python environment -------------
# Pin RETICULATE_PYTHON directly to the env's interpreter (must be set before reticulate initializes any Python session).
# This is more robust than use_condaenv(name) on this machine, since that relies on `conda env list` to discover environments, which does not reliably see micromamba envs.
# Adjust the path if your pan-human-azimuth env lives elsewhere (check with: micromamba env list).
Sys.setenv(RETICULATE_PYTHON = "./.conda-pan-human-azimuth/bin/python")

# Confirm reticulate actually loaded THIS interpreter and that it can see panhumanpy inside it. 
py_config()
# might expose token
# Sys.getenv()
Sys.getenv(c("RETICULATE_PYTHON", "CONDA_PREFIX", "CONDA_DEFAULT_ENV"))

cat("Active Python:", py_config()$python, "\n")
cat("panhumanpy version:", as.character(import("panhumanpy")$`__version__`), "\n")

# 3. Run Pan-human Azimuth annotation locally --------------------------------
# output_mode = 'detailed' additionally returns the raw per-level labels (level_1_labels ... level_N_labels), used in 02_hierarchy_consistency_check.R to inspect exactly where a cell's hierarchy becomes inconsistent. First run downloads and caches the model weights locally - can take a while.
bmcite_ann <- ANNotate(bmcite, output_mode = "detailed")
saveRDS(bmcite_ann, file.path(out_obj_dir, "bmcite_ann.rds"))

# 4. Examine model output -----------------------------------------------------
md <- bmcite_ann@meta.data
colnames(md)

# save the metadata
write_tsv(md,file = "../out/table/bmcite_ann.tsv")

# Key columns:
# - full_hierarchical_labels : combined "A|B|C|..." path for each cell
# - full_consistent_hierarchy: TRUE/FALSE, whether that path is internally consistent across all predicted levels
# - final_level_labels       : label at the deepest level reached (depth varies per cell - not a fixed granularity)
# - final_level_confidence   : softmax/calibrated confidence at that final level
# - azimuth_broad/medium/fine: post-processed labels at fixed granularity; "False" for any cell with full_consistent_hierarchy == FALSE

table(md$full_consistent_hierarchy)

ggplot(md, aes(x = final_level_confidence)) +
  geom_histogram(bins = 20, fill = "skyblue") +
  labs(x = "Softmax probability", y = "Count",
       title = "Histogram of final-level confidence") +
  theme_bw()

# Confidence vs. hierarchy-consistency (as discussed in the upstream issue)
ggplot(md, aes(x = full_consistent_hierarchy, y = final_level_confidence)) +
  geom_boxplot() +
  labs(x = "full_consistent_hierarchy", y = "final_level_confidence") +
  theme_bw()

# QC / filtering step, as suggested in the AzimuthAPI vignette
bmcite_qc <- subset(bmcite_ann, final_level_confidence > 0.5)

# 5. UMAP of Azimuth embeddings -----------------------------------------------
bmcite_qc <- RunUMAP(bmcite_qc, dims = 1:128, reduction = "azimuth_embed",
                      reduction.name = "azimuth_umap")

p_full <- DimPlot(bmcite_qc, group.by = "full_hierarchical_labels",
                   label = TRUE, label.size = 1.5, reduction = "azimuth_umap",
                   repel = TRUE) + NoLegend()
p_final <- DimPlot(bmcite_qc, group.by = "final_level_labels",
                    label = TRUE, label.size = 2, reduction = "azimuth_umap") +
  NoLegend()
p_conf <- FeaturePlot(bmcite_qc, features = "final_level_confidence",
                       reduction = "azimuth_umap")

ggsave(file.path(out_plot_dir, "azimuth_full_hierarchical_labels.png"), p_full,  width = 7, height = 6)
ggsave(file.path(out_plot_dir, "azimuth_final_level_labels.png"),      p_final, width = 7, height = 6)
ggsave(file.path(out_plot_dir, "azimuth_final_level_confidence.png"),  p_conf,  width = 7, height = 6)

# 6. Refined labels at consistent granularity (added automatically by ANNotate() by default: azimuth_broad / azimuth_medium / azimuth_fine) --
p_medium <- DimPlot(bmcite_qc, group.by = "azimuth_medium",
                     label = TRUE, label.size = 3, reduction = "azimuth_umap") +
  NoLegend()
p_fine <- DimPlot(bmcite_qc, group.by = "azimuth_fine",
                   label = TRUE, label.size = 3, reduction = "azimuth_umap") +
  NoLegend()

ggsave(file.path(out_plot_dir, "azimuth_medium.png"), p_medium, width = 7, height = 6)
ggsave(file.path(out_plot_dir, "azimuth_fine.png"),   p_fine,   width = 7, height = 6)

# Collapse rare fine labels (< 20 cells) for a cleaner visualization
bmcite_qc <- PrepLabel(bmcite_qc, "azimuth_fine", "azimuth_fine_filtered", cutoff = 20)
p_fine_filtered <- DimPlot(bmcite_qc, group.by = "azimuth_fine_filtered",
                            label = TRUE, label.size = 3, reduction = "azimuth_umap") +
  NoLegend()
ggsave(file.path(out_plot_dir, "azimuth_fine_filtered.png"), p_fine_filtered, width = 7, height = 6)

saveRDS(bmcite_qc, file.path(out_obj_dir, "bmcite_qc.rds"))

# 7. (Optional) QC expression heatmaps by predicted cell type ----------------
plots <- make_azimuth_QC_heatmaps(bmcite_qc)
length(plots)

plots[["Immune_Lymphoid cell_1"]]
plots[["Immune_Myeloid cell_1"]]
