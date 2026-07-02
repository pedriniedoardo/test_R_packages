## ----installs, eval=FALSE-----------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
# # orthogene is only available on Bioconductor>=3.14
# if(BiocManager::version()<"3.14")
#   BiocManager::install(update = TRUE, ask = FALSE)
# 
# BiocManager::install("orthogene")

## ----setup--------------------------------------------------------------------
library(orthogene)

data("exp_mouse")
# Setting to "homologene" for the purposes of quick demonstration.
# We generally recommend using method="gprofiler" (default).
method <- "homologene"  

## ----infer_species------------------------------------------------------------
# use a matrix of expression
matches <- orthogene::infer_species(gene_df = exp_mouse, 
                                    method = method)

# use a vector
GOI <- rownames(exp_mouse)
orthogene::infer_species(GOI,method = method)

## ----convert_orthologs--------------------------------------------------------
exp_rat <- orthogene::convert_orthologs(gene_df = exp_mouse, 
                                        input_species = "mouse", 
                                        output_species = "rat",
                                        method = method)

## ----infer_species2-----------------------------------------------------------
matches <- orthogene::infer_species(gene_df = exp_rat, 
                                    method = method)

# use a vector
GOI <- rownames(exp_rat)
orthogene::infer_species(GOI,method = method)

## ----convert_orthologs2-------------------------------------------------------
exp_human <- orthogene::convert_orthologs(gene_df = exp_mouse, 
                                          input_species = "mouse", 
                                          output_species = "human",
                                          method = method)

## ----infer_species3-----------------------------------------------------------
matches <- orthogene::infer_species(gene_df = exp_human, 
                                    method = method)

# use a vector
GOI <- rownames(exp_human)
orthogene::infer_species(GOI,method = method)

## ----infer_species4-----------------------------------------------------------
matches <- orthogene::infer_species(gene_df = exp_human, 
                                    test_species = method, 
                                    method = method)

## ----Session Info-------------------------------------------------------------
utils::sessionInfo()
