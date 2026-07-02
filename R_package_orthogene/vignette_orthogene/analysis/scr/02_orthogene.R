## ----installs, eval=FALSE-----------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
# # orthogene is only available on Bioconductor>=3.14
# if(BiocManager::version()<"3.14") BiocManager::install(version = "3.14")
# 
# BiocManager::install("orthogene")

## ----setup--------------------------------------------------------------------
library(orthogene)

data("exp_mouse")
# Setting to "homologene" for the purposes of quick demonstration.
# We generally recommend using method="gprofiler" (default).
method <- "homologene"  

## ----convert_orthologs--------------------------------------------------------
gene_df <- orthogene::convert_orthologs(gene_df = exp_mouse,
                                        gene_input = "rownames", 
                                        gene_output = "rownames", 
                                        input_species = "mouse",
                                        output_species = "human",
                                        non121_strategy = "drop_both_species",
                                        method = method) 

knitr::kable(as.matrix(head(gene_df)))

## ----map_species--------------------------------------------------------------
species <- orthogene::map_species(species = c("human",9544,"mus musculus",
                                              "fruit fly","Celegans"), 
                                  output_format = "scientific_name")
print(species)

## ----report_orthologs---------------------------------------------------------
orth_zeb <- orthogene::report_orthologs(target_species = "zebrafish",
                                        reference_species = "human",
                                        method_all_genes = method,
                                        method_convert_orthologs = method) 
knitr::kable(head(orth_zeb$map))
knitr::kable(orth_zeb$report)

## ----map_genes----------------------------------------------------------------
genes <-  c("Klf4", "Sox2", "TSPAN12","NM_173007","Q8BKT6",9999,
            "ENSMUSG00000012396","ENSMUSG00000074637")
mapped_genes <- orthogene::map_genes(genes = genes,
                                     species = "mouse", 
                                     drop_na = FALSE)
knitr::kable(head(mapped_genes))

## ----aggregate_mapped_genes---------------------------------------------------
data("exp_mouse_enst") 
knitr::kable(tail(as.matrix(exp_mouse_enst)))

exp_agg <- orthogene::aggregate_mapped_genes(gene_df=exp_mouse_enst,
                                             input_species="mouse", 
                                             agg_fun = "sum")
knitr::kable(tail(as.matrix(exp_agg)))

## ----all_genes----------------------------------------------------------------
genome_mouse <- orthogene::all_genes(species = "mouse", 
                                     method = method)

knitr::kable(head(genome_mouse))

## ----Session Info-------------------------------------------------------------
utils::sessionInfo()
