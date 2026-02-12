# Creation of a new folder in the repository<br>
# open your repository, go to your folder of interest
# - create a new folder<br>
# - copy and paste in this new folder: 
  
# - metadata.csv created previously: to update it, open metadata.csv in Excel, fill in the file and save it using .csv with “;” as separator.
# - rRMSAnalyzer_template.R (https://github.com/RibosomeCRCL/rRMSAnalyzer/tree/main/vignettes)
# - "human.methylated.rda" and "human.suspected.rda"  (https://github.com/RibosomeCRCL/rRMSAnalyzer/tree/main/data)
# - a folder with all the 5’end read counts files
  
# Perform analysis using RStudio
## Open RStudio
## Update RMSAnalyzer
library(rRMSAnalyzer)

# Create a RiboClass (metadata, counts, C-score parameters)
# As an example, the provided dataset "ribo_toy" is used.
data("ribo_toy")

# Personal data preparation
## Prepare the comparison table
# - open the metadata.csv using Excel
# - add as many new columns as required (1 column per comparison, here "column name" = "variableA to C")
# - complete the new columns with
# - ctrl = control condition
# - P1 = condition 1
# - P2 = condition 2….
# • NA = not applicable
# Save your metadata file using .csv and separator ";"
# In Rstudio, load a newRiboClass using:<br>

# ribo <- load_ribodata(
#   #data & metadata files path
#   count_path = "~/Desktop/RMS/path/counts/",
#   metadata = "~/Desktop/RMS/path/metadata.csv",# data & metadata files separator
#   count_sep = "\t",
#   metadata_sep = ";",
#   # count data parameters :
#   count_header = FALSE,
#   count_value = 3,
#   count_rnaid = 1,
#   count_pos = 2,
#   # Metadata parameters :
#   metadata_key = "filename",
#   metadata_id = "samplename",
#   # c-score parameters :
#   flanking = 6,
#   method = "median",
#   ncores = 1)


ribo <- ribo_toy


# Check importation by ploting a PCA and COA ------------------------------
# change column name of the lib if necessary #ribo_toy lib_col = "run"
plot_pca(ribo,"run")

plot_coa(
  ribo,
  color_col = NULL,
  axes = c(1, 2),
  only_annotated = FALSE,
  title = "default",
  subtitle = "default",
  object_only = FALSE
)


# Perfom QC analysis ------------------------------------------------------
### Generate a QCreport
# ribo_toy lib_col = "run"
report_qc(ribo = ribo, specie = "human", library_col = "run", project_name = "ribo_toy")

### Adjustment using ComBat-Seq if necessary
# from the c score pca it is quite clear the batch effect of the run
ribo_adj <- adjust_bias(ribo_toy,"run") 

# verify the adjustment using a PCA
# ribo_toy lib_col = "run"
plot_pca(ribo_adj,"run")

# ribo_toy lib_col = "run"
report_qc(ribo = ribo_adj, specie = "human", library_col = "run", project_name = "ribo_toy_adj", comments = "./comment_QC.Rmd")


# Remove problematic or unused samples (optional) -------------------------
# To keep a small number of samples
# ribo_adj_small <- keep_ribo_samples(ribo_adj,c("sample1","sample2","sample3","sample4"))

# To remove a small number of samples
ribo_adj_small <- remove_ribo_samples(ribo_adj,c("S7","RNA1","RNA2"))

### Uniformisation of the RiboClass name
ribo_adj <- ribo_adj_small


# Annotation of human rRNA 2’Ome sites ------------------------------------
# Choose between 2 annotations: known (n=112 sites) or suspected (n=17). This section or the next one. If the two annotations are required, both of known and suspected sites analysis are desired, the 2 scripts has to be run twice from here.

# Annotation of the 112 known sites 
data("human_methylated")
cat("human_methylated's rna names: ", unique(human_methylated$rRNA),"\n")
cat("ribo's rna names: ", as.character(ribo$rna_names$current_name)) # change name of RiboClass if necessary

# change name of RiboClass if necessary for ribo_adj if you adjusted the data
ribo_adj_name <- rename_rna(ribo_adj,
                            new_names = c("5S", "5.8S", "18S", "28S"))

ribo_adj_annot <- annotate_site(ribo_adj_name,
                                annot = human_methylated,
                                anno_rna = "rRNA",
                                anno_pos = "Position",
                                anno_value = "Nomenclature")

# Annotation of the 19 suspected sites
# load(file = "/home/path/to/RMS/human.suspected.rda") 
data("human_suspected")

# verify rRNA names
ribo_adj_annot <- annotate_site(ribo_adj_annot, human_suspected)

ribo_adj_annot <- annotate_site(ribo_adj_annot,
                                annot = human_suspected,
                                anno_rna = "rRNA",
                                anno_pos = "Position",
                                anno_value = "Nomenclature")


# Analytic reports --------------------------------------------------------
### Reports for comparisons explained in column 1
# In the previous example: in column “VariableA”, comparison ctrl vs P1 and ctrl vs P2
# notuce that I need to generate in the metadata as many columns as comparisons. fill with NA the entries that do not have a match for the comparison

#### Filter the data
# To do only if you want to compare a subset of your samples (i.e., not all the samples included in the RiboClass).
kept_samples <- ribo_adj_annot$metadata %>%
  # keep lines that are not "NA"
  # dplyr::filter(!is.na(variableA)) %>% #ribo_toy (!is.na(comp1))
  dplyr::filter(!is.na(comp1)) %>%
  dplyr::pull(samplename)

# create a new RiboClass including the subdata variableA
# ribo_toy ribo_adj_annot_comp1 <- ...
ribo_adj_annot_variableA <- keep_ribo_samples(ribo_adj_annot,kept_samples) 

# Create mandatory comparison table for the report_diff_sites  ------------
comparisons <- tibble::tibble(
  # name of the comparisons, add as many as required
  # comp = c("comp1", "comp2"), 
  comp = c("comp1"), 
  # name of the ctrl as appearing in your metadata for each comparison, in the order, add as many name than comparison that will be performed
  # ctrl = c("ctrl", "ctrl"),
  ctrl = c("cond1"),
  # name of the cases as appearing in your metadata for each comparison, in the order, add as many name than comparison that will be performed
  # cases = c("P1", "P2")
  cases = c("cond2")
)


# Generate a report to compare 2’Ome profile    ---------------------------
# condition_col: change the name of the column to perform new analysis. ribo_toy condition_col = "comp1"
# comments: change the path and name of your files with your comments if required
report_2ome_sites(ribo = ribo_adj_annot_variableA,
                  specie = "human",
                  # condition_col = "variableA",
                  condition_col = "comp1",
                  project_name = "test",
                  comments = "../out/comment/comment_2ome_variableA.Rmd") 

# Generate a report to perform site-by-site 2’Ome comparison --------------
# ribo_toy : ribo = ribo_adj_annot_comp1
# condition_col: change the name of the column to perform new analysis . ribo_toy condition_col = "comp1"
# comparisons: name of the comparison table
# comments: change the name of your files with your comments if required
report_diff_sites(ribo = ribo_adj_annot_variableA,
                  specie = "human",
                  # condition_col = "variableA",
                  condition_col = "comp1",
                  project_name = "test",
                  comparisons = comparisons,
                  comments = "../out/comment/comment_diff_site_variableA.Rmd")

# Reports for comparisons explained in column 2 ---------------------------
# Re-do all the steps of [Reports for comparisons explained in column 1] but using the name of your column of interest


# Export data -------------------------------------------------------------
# change for the annotated RiboClass if needed
ribo_df <- extract_data(ribo_adj_annot,
                        col = "cscore",
                        only_annotated = TRUE)

# give path
write.csv(ribo_df, "../out/table/cscore_df.csv")

# Check local environments of sites of interest (optionnal)  --------------
# If you want to verify a site that shows up as significant in the last report
# For all samples

plot_counts_env(ribo_adj_annot, rna= "18S", pos = 1442, samples = "all", condition = "condition")

# For some samples
plot_counts_env(ribo_adj_annot, rna = "18S", pos = 1442, samples = c("S1", "S2"))
