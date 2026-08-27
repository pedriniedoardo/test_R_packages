# AzimuthAPI - package installation reference -------------------------------
# This project uses renv. Install with renv (not install.packages/remotes
# directly) so the packages land in the project-local renv library and get
# captured by renv::snapshot().
#
# AzimuthAPI itself is not on CRAN/Bioconductor - it and SeuratData are
# installed straight from GitHub.

# renv::install("remotes")
# renv::install("reticulate")
# renv::install("tidyverse")
# renv::install("satijalab/seurat-data")
# renv::install("Seurat")
# renv::install("satijalab/AzimuthAPI")

# After installing, snapshot the lockfile:
# renv::snapshot()

# Local (non-renv) requirement: a Python environment with the `panhumanpy`
# package, used via reticulate for the ANNotate() local-inference path
# (see 01_test_AzimuthAPI_local.R). Set up once, outside of R, e.g.:
#
#   conda create -n pan-human-azimuth python=3.9
#   conda activate pan-human-azimuth
#   pip install panhumanpy
#
# renv does not manage this Python environment; reticulate just points at it.
