# Working BisqueRNA Example
# This code walks through generating sample data and then applying Bisque's two main decomposition functions.


# 1. Setup ----------------------------------------------------------------
# First, you'll need to install BisqueRNA and its dependency Biobase from Bioconductor.

# Load the required libraries
library(BisqueRNA)
library(Biobase)

# Set a seed for reproducibility
set.seed(42)

# 2. Simulate Data --------------------------------------------------------
# We will create a dataset with both single-cell and bulk RNA-seq data for 10 individuals. The SimulateData function is a helper from the package that makes this easy. The output includes the single-cell ExpressionSet (sc.eset), the bulk ExpressionSet (bulk.eset), the true cell type proportions (props), and a list of marker genes (markers).

# To mimic a real-world scenario, we will then remove the single-cell data for individuals 1-5, leaving only their bulk data to be decomposed.

# Define cell types and their average proportions
cell.types <- c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelial Cells")
avg.props <- c(0.5, 0.2, 0.2, 0.07, 0.03)

# Generate the full dataset for 10 individuals
sim.data <- SimulateData(n.ind = 10, n.genes = 100, n.cells = 500, 
                         cell.types = cell.types, avg.props = avg.props)

# Subset the single-cell data to only include individuals 6-10
# This creates a reference set that does not contain all bulk samples.
sc.eset <- sim.data$sc.eset[, sim.data$sc.eset$SubjectName %in% as.character(6:10)]

# The bulk data still contains all 10 individuals
bulk.eset <- sim.data$bulk.eset

# Store the ground truth proportions and marker genes for later validation
true.props <- sim.data$props
markers <- sim.data$markers

# Clean up the large simulation object
rm(sim.data)

print("Data simulation complete.")


# 3. Reference-Based Decomposition ----------------------------------------
# This is the primary mode for Bisque. It uses the complete expression profiles from a single-cell reference (sc.eset) to estimate cell type proportions in the bulk RNA-seq samples (bulk.eset). We set use.overlap=TRUE because our reference (individuals 6-10) contains some of the same individuals present in our bulk dataset.

# Run the reference-based decomposition
res.ref <- BisqueRNA::ReferenceBasedDecomposition(
  bulk.eset = bulk.eset,
  sc.eset = sc.eset,
  markers = NULL,       # Using all genes, not just markers
  use.overlap = TRUE    # Specify that some individuals have both scRNA and bulk data
)

# The estimated proportions are in the 'bulk.props' slot
ref.based.estimates <- res.ref$bulk.props

# Print the estimated proportions
# cat("## Reference-Based Decomposition Results:\n")

print(round(ref.based.estimates, 3))

# Validate the results by correlating estimates with the true proportions
r <- cor(as.vector(ref.based.estimates), 
         as.vector(true.props[rownames(ref.based.estimates), colnames(ref.based.estimates)]))

# 4. Marker-Based Decomposition -------------------------------------------
# This mode is useful when you don't have a full single-cell reference dataset but know which genes are markers for your cell types. It uses the expression of only these marker genes to infer the relative abundance of each cell type.

# Important: The resulting proportions are relative within each cell type (i.e., across samples), not absolute. You cannot directly compare the abundance of Neurons to Astrocytes in a given sample.

# Run marker-based decomposition using the pre-defined markers
res.marker <- BisqueRNA::MarkerBasedDecomposition(
  bulk.eset = bulk.eset, 
  markers = markers, 
  weighted = FALSE
)

# The estimated proportions are in the 'bulk.props' slot
marker.based.estimates <- res.marker$bulk.props

# Print the estimated relative abundances

print(round(marker.based.estimates, 3))

# To validate, we must scale the true proportions for each cell type
# This makes them relative, just like the marker-based estimates
scaled.true.props <- t(scale(t(true.props)))[rownames(marker.based.estimates),]

# Correlate the marker-based estimates with the scaled true proportions
r.marker <- cor(as.vector(marker.based.estimates),
                as.vector(scaled.true.props))
