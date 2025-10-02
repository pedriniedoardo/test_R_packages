# AIM ---------------------------------------------------------------------
# test the tool with a toy dataset to make a fast and reproducible example

# libraries ---------------------------------------------------------------
library(Seurat)
library(BayesPrism)
library(tidyverse)

# Step 1: Setup and Data Simulation ---------------------------------------
# create some mock data. We'll simulate a Seurat object to act as our single-cell reference and a simple matrix for our bulk RNA-seq samples.

# Simulate a Structured scRNA-seq Reference
set.seed(123)
n_genes <- 2000
n_cells <- 2000 # Using more cells for a better reference
cell_types <- c("T-cell", "B-cell", "Macrophage", "Endothelial")

# Assign cells to types
cell_type_labels <- sample(cell_types, n_cells, replace = TRUE)
names(cell_type_labels) <- paste0("Cell_", 1:n_cells)
gene_names <- paste0("Gene", 1:n_genes)

# Define marker genes for each cell type
marker_gene_list <- list(
  `T-cell`      = gene_names[1:50],
  `B-cell`      = gene_names[51:100],
  `Macrophage`  = gene_names[101:150],
  `Endothelial` = gene_names[151:200]
)

# Create a base matrix of mean expression values (lambda)
# Low background expression for all genes
lambda_matrix <- matrix(0.1, nrow = n_genes, ncol = length(cell_types),
                        dimnames = list(gene_names, cell_types))

# Increase expression for marker genes within their respective cell type
for (cell_type in cell_types) {
  lambda_matrix[marker_gene_list[[cell_type]], cell_type] <- 8.0
}

# Generate the final count matrix using these mean expression values
sc_counts <- sapply(1:n_cells, function(i) {
  current_cell_type <- cell_type_labels[i]
  mean_expression_vector <- lambda_matrix[, current_cell_type]
  
  # Use rpois to generate counts based on the defined expression pattern
  counts <- rpois(n_genes, lambda = mean_expression_vector)
  return(counts)
})

# rpois(n_genes,lambda = lambda_matrix[, "Endothelial"])
# dimnames(sc_counts) <- list(gene_names, names(cell_type_labels))

# Create the Seurat object
seurat_ref <- CreateSeuratObject(counts = sc_counts)
seurat_ref$cell_type <- cell_type_labels

# save the sc object
saveRDS(seurat_ref, file = "../out/object/seurat_ref.rds")

# produce the aggregared counts per cell type
cts <- AggregateExpression(object = seurat_ref,
                           group.by = c("cell_type"),
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = FALSE)

saveRDS(cts$RNA %>% as.matrix(),file = "../out/object/seurat_ref_pBulkCellType.rds")

# Simulate a pbluk by changing the proportion of different cell types

# split the dataset into pure cell types to use for sampling
Idents(seurat_ref) <- "cell_type"

list_sc_pure <- lapply(cell_types,function(x){
  subset(seurat_ref,subset = cell_type == x)
}) %>%
  setNames(cell_types)

# now sample from different proportions of cells load the lut
# to avoid issues in the selection of the cells, add two cells to all the categoreis
list_LUT_prop <- read_csv("../data/LUT_prop_pBulk.csv") %>%
  pivot_longer(names_to = "cell_type",values_to = "prop",-sample_id) %>%
  mutate(n_cells_1p = round(450*prop)+2) %>%
  split(.,f = .$sample_id)

# split the data per pBulk sample
# samp <- list_LUT_prop$Sample_1_T_dom

# celltype <- samp$cell_type[1]
# ncell <- samp$n_cells_1p[1]
# ncell <- 2


set.seed(123)

list_pBulk_test <- list_LUT_prop %>%
  lapply(., function(samp){
    # check the progress
    print(samp)
    
    # pull the name of the sample
    samp_name <- unique(samp$sample_id)
    
    # select the cells to generate the pBulks
    pBulk_sample <- pmap(list(samp$cell_type,samp$n_cells_1p),function(celltype,ncell){
      # select the specific pool of cells from the pure datasets
      sc_test <- list_sc_pure[[celltype]]
      
      # select at random the number of cells
      random_cells <- sample(size = ncell,x = rownames(sc_test@meta.data),replace = F)
      sc_test$test <- rownames(sc_test@meta.data) %in% random_cells
      
      # confirm the number of cells
      sum(sc_test$test)
      
      # subset the selected cells
      sc_test_final <- subset(sc_test,subset = test == 1)
      
      # generate the pBulk
      cts <- AggregateExpression(object = sc_test_final,
                                 group.by = c("cell_type"),
                                 assays = 'RNA',
                                 slot = "counts",
                                 return.seurat = FALSE)
      
      data <- cts$RNA %>% as.data.frame()
      return(data)
    }) %>%
      bind_cols()
    
    # aggregate the reads in one sample
    df_sample <- rowSums(pBulk_sample) %>%
      as.data.frame(.) %>%
      rownames_to_column("gene")
    
    colnames(df_sample)[2] <- samp_name
    
    return(df_sample)
  })

# make a unique table with all the samples and convert it to a matrix
bulk_counts <- list_pBulk_test %>%
  reduce(full_join,by=c("gene")) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# save the siumulated matrix of pBulks samples
saveRDS(bulk_counts, file = "../out/object/bulk_counts.rds")

# Check the dimensions of our simulated data
print(dim(seurat_ref))
print(dim(bulk_counts))

#  Step 2: Prepare Data for BayesPrism ------------------------------------
# BayesPrism requires the data in specific formats. We need to extract the raw counts and cell labels from our Seurat object and ensure the gene names match between our reference and bulk data.

# Extract Information from the Seurat Object
# BayesPrism works best with raw counts.
sc_counts_matrix <- GetAssayData(seurat_ref, slot = "counts")

# Get cell type labels. These will be our 'cell.type.labels'.
cell_types <- seurat_ref$cell_type

# Get cell state labels. For a simple case, these can be the same as cell types.
# In more complex analyses, these could be cluster IDs or subtypes.
# https://github.com/Danko-Lab/BayesPrism/issues/66
cell_states <- seurat_ref$cell_type

# Filter and Align Genes
# reshape the data to be accepted by the tool
# generally data are reporeted as feature X sample, but in the tool they are handled as sample X feature.
# we can input sparse matrices: https://github.com/Danko-Lab/BayesPrism/issues/58
sc.dat <- t(sc_counts_matrix)
sc.dat[1:5,1:5]
dim(sc.dat)

bk.dat <- t(bulk_counts)
bk.dat[1:5,1:5]
dim(bk.dat)

# Run the reccommended filtering on the genes in this step for a real data analysis

# Step 3: Run BayesPrism Deconvolution ------------------------------------
# With the data properly formatted and aligned, we can now run the deconvolution.

# Construct the BayesPrism Object
# We combine the reference and cell labels into a 'prism' object.
my_prism <- new.prism(
  reference = sc.dat,
  mixture = bk.dat,
  key = NULL, # No specific tumor/normal key needed for this example
  cell.type.labels = seurat_ref$cell_type,
  cell.state.labels = seurat_ref$cell_type,
  input.type = "count.matrix"
)

# Run the Deconvolution Algorithm
# This is the core computational step.
# We'll use 4 cores for this example. Adjust 'n.cores' as needed.
bp_results <- run.prism(prism = my_prism, n.cores = 20)

# You can look at the returned object to see what it contains
bp_results

# Step 4: Extract and Calculate Deconvoluted Reads ------------------------
# this can be extracted from the object as follows
# x <- "Endothelial"

list_bulk <- lapply(unique(cell_types), function(x){
  mat <- get.exp(bp=bp_results,
                 state.or.type="type",
                 cell.name=x)
  
  data <- mat %>%
    t()
  
  return(data)
}) %>%
  setNames(unique(cell_types))

# generate integer estimates
list_bulk_integer <- lapply(unique(cell_types), function(x){
  mat <- get.exp(bp=bp_results,
                 state.or.type="type",
                 cell.name=x)
  
  # generate integer level estimates
  # https://github.com/Danko-Lab/BayesPrism/issues/121
  # round up the Z matrix to the nearest integer and use it as input for DESeq2 (without any normalization step).
  data <- mat %>%
    t() %>%
    round()
  return(data)
}) %>%
  setNames(unique(cell_types))

# get the fraction estimated
theta <- get.fraction(bp=bp_results,
                      which.theta="final",
                      state.or.type="type")
theta

# test --------------------------------------------------------------------
# compare the expected prop vs the actula ones
df_cor <- full_join(
  read_csv("../data/LUT_prop_pBulk.csv") %>%
    pivot_longer(names_to = "cell_type",values_to = "prop",-sample_id),
  theta %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    pivot_longer(names_to = "cell_type",values_to = "prop",-sample_id),by = c("sample_id","cell_type"),suffix = c(".expexted",".BayesPrism"))

df_cor %>%
  ggplot(aes(x=prop.expexted,y=prop.BayesPrism)) +
  geom_point() +
  geom_abline(intercept = 0,slope = 1,col="red",linetype = "dashed") +
  theme_bw() +
  facet_wrap(~cell_type)+
  theme(strip.background = element_blank())

# confirm that the sum of the deconvoluted pBulk per sample correponds to the actual ionput values
# bulk_convoluted <- list_bulk$Endothelial + list_bulk$`T-cell` + list_bulk$`B-cell`+ list_bulk$Macrophage
bulk_convoluted <- purrr::reduce(list_bulk,`+`)

# they are the same
all.equal(bulk_convoluted,bulk_counts)
