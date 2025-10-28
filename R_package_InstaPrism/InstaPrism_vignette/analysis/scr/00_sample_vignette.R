# Introduction ------------------------------------------------------------
# InstaPrism is an R package for cell type composition and gene expression deconvolution in bulk RNA-Seq data.

# Based on the same conceptual Bayesian framework as BayesPrism, InstaPrism re-implements BayesPrism in a derandomized framework by replacing the time-consuming Gibbs sampling steps in BayesPrism with a fixed-point algorithm, which greatly accelerated the calculation speed while maintaining highly comparable performance. It works as an independent R package and does not require the users to have BayesPrism installed.

# This tutorial is associated the latest version of InstaPrism (InstaPrism v0.1.6). extdata folder associated with this tutorial can be downloaded from zenodo repository.

# Quick start with your own data
# Load the InstaPrism package.

library(InstaPrism)
library(ggpmisc)

# 1. Prepare the bulk input -----------------------------------------------
# Here we provide an example of TPM-normalized bulk expression data from ovarian cancer. For deconvolution with InstaPrism, the bulk input should be in a non-log transformed scale (TPM or count).
bulk_expr <- read.csv(system.file('extdata','example_bulk.csv',package = 'InstaPrism'))
saveRDS(bulk_expr,file = "../data/bulk_expr.rds")
bulk_expr[1:5,1:5]

# 2. Prepare the reference ------------------------------------------------
# A reference profile containing prior knowledge from the same tissue type is required for InstaPrism deconvolution. With the provided ovarian cancer bulk example, we can utilize the built-in ovarian cancer reference.

OV_ref <- InstaPrism_reference('OV') 
# This reference contains two layers of information: 1) a probability matrix phi.cs that specifies the event probability of a gene being expressed in a cell state; and 2) a map list that specifies the corresponding cell type for each cell state.


# 3. Deconvolution with InstaPrism ----------------------------------------
# Upon preparation of the necessary inputs, the deconvolution task can be accomplished with the InstaPrism() function:
deconv_res <- InstaPrism(bulk_Expr = bulk_expr, refPhi_cs = OV_ref)

# 4. Obtain deconvolution results -----------------------------------------
# InstaPrism provides both cell type fraction and gene expression deconvolution. The deconvolved cell type fraction is accessible with:
estimated_frac <- t(deconv_res@Post.ini.ct@theta)
head(estimated_frac)

# The cell types in the deconvolved cell type fraction table correspond to the cell types provided in the reference. In other words, given the reference, InstaPrism will return posterior estimates for the cell types in the reference.

# The cell-type specific gene expression deconvolution results can be obtained by
Z <- get_Z_array(deconv_res)

# This Z array is a sample by gene by cell-type array
dim(Z)

# To extract a cell-type specific gene expression, one can use

mal_Z <- Z[,,'malignant']
head(mal_Z[,1:5])

# Or directly use reconstruct_Z_ct_initial() function to get a specific gene expression of interest
mal_Z <- reconstruct_Z_ct_initial(InstaPrism_obj = deconv_res, cell.type.of.interest = 'malignant')

# 5. Examine deconvolution performance ------------------------------------
# We provide a useful deconv_performance_plot() function to visualize the correlation between the estimated and ground truth cell type fractions. Since the names of the estimated and ground truth cell types may not always match, the function incorporates an automatic matching step. Specifically, for each cell type in the ground truth, the function identifies the most strongly correlated cell type in the estimated fractions, and labeled it as maxCorName in the visualization.

# load ground truth cell type fractions for the example bulk data
truth_frac <- read.csv(system.file('extdata','truth_frac.csv',package = 'InstaPrism'))
deconv_performance_plot(truth_frac = truth_frac,
                        estimated_frac = t(deconv_res@Post.ini.ct@theta))



# De novo reference construction ------------------------------------------
# In the previous example we demonstrated a deconvolution process using the InstaPrism built-in reference. In real practice, users can generate their own reference profiles using the refPrepare() function, which takes the following inputs:

# scExpr: a single-cell expression data as prior information (raw count or CPM normalized expression)
# cell.type.labels: a character vector indicating cell types of each cell from the scRNA data
# cell.state.labels: a character vector indicating cell states of each cell from the scRNA data

# Cell type labels indicate the broad classification of cells based on their biological or functional identity, such as T cells, B cells, or malignant cells; whereas cell state labels indicate more specific statuses within each cell type, which can be derived from subclustering information from each cell type. Accounting for cell state specifications is crucial for addressing the heterogeneity within a cell type, thereby enhancing the robustness of deconvolution results.

# Below, we provide an example showing how to construct a reference with user-provided single-cell data.


# 1. Toy example ----------------------------------------------------------
data("sim.data")
# sim.data is a list generated using the BisqueRNA::SimulateData() function, it contains the following objects:

# sc.eset: A simulated single cell expression object in raw count format, with 1000 genes and 10,000 single cells, the cells are annotated with sample ID and cell type labels
# bulk.eset: A simulated bulk expression object for 20 samples
# prop: A ground truth cell-type fraction table for 20 simulated bulk samples
# markers: A dataframe indicating marker genes for each cell type

library(Biobase)
sc.eset <- sim.data$sc.eset
scExpr <- exprs(sc.eset)
cell_type_labels <- pData(sc.eset)$cellType

# In this toy example, since the single-cell data is simulated to reflect only differences at the cell-type level, without additional variations within the cell types during the simulation, we can skip cell state identification and simplifies the reference construction as follows:
refPhi_obj <- refPrepare(sc_Expr = scExpr,
                         cell.type.labels = cell_type_labels,
                         cell.state.labels = cell_type_labels)

# This reference object can then be passed to InstaPrism() function for deconvolution
toy_example_res <- InstaPrism(bulk_Expr = exprs(sim.data$bulk.eset),refPhi_cs = refPhi_obj)
deconv_performance_plot(truth_frac = t(sim.data$props),estimated_frac = t(toy_example_res@Post.ini.ct@theta))

# We also provide a sample code demonstrating how simdata is created.
library(Biobase)
library(BisqueRNA)

cell.types <- c("Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Endothelial Cells","others")
avg.props <- c(.5, .2, .15, .07, .03,.05)
sim.data <- SimulateData(n.ind=20,
                         n.genes=1000,
                         n.cells=500,
                         cell.types=cell.types,
                         avg.props=avg.props)

# 2. Cell state specification ---------------------------------------------
# In real practice, constructing a reference from single-cell data is more complex, and one cell type may not be represented by a single profile. Ideally, if finer subclustering information is available, it can directly serve as cell state labels. However, when such information is missing, we provide a get_subcluster() function to obtain detailed cell state information within each cell type. This function supports two types of subclustering methods: scran and SCISSORS.

# scran
# This method implements the quickCluster() function from scran package for fast identification of subclustering information. It includes a min.subcluster.size hyperparameter to specify the minimum size of each subcluster.
library(scran)

cell_state_lables_scran <- get_subcluster(scExpr = scExpr,
                                          cell_type_labels = cell_type_labels,
                                          subcluster_method = 'scran',
                                          min.subcluster.size = 200)
table(cell_state_lables_scran)

# SCISSORS
# This method implements the SCISSORS method for sub-clusters identification. It contains a SilhouetteScores.cutoff hyperparamter to determine which cell types to be re-culstered.
library(biomaRt)
library(phateR)
library(SCISSORS)

cell_state_lables_SCISSORS <- get_subcluster(scExpr = scExpr,
                                             cell_type_labels = cell_type_labels,
                                             subcluster_method = 'SCISSORS',
                                             SilhouetteScores.cutoff = 0.7)
table(cell_state_lables_SCISSORS)

# Note that for this example data, using the current SilhouetteScores.cutoff, we received the message: “No cell types need to be reclustered under the given SilhouetteScores.cutoff; original cell type labels will be returned.” This indicates that the cell types already demonstrate good internal consistency and are distinct from each other, thus no further subclustering is needed.

# The identified cell state labels can be directly pass to refPrepare() function:
refPhi_obj1 <- refPrepare(sc_Expr = scExpr, cell.type.labels = cell_type_labels, cell.state.labels = cell_state_lables_scran)
refPhi_obj2 <- refPrepare(sc_Expr = scExpr, cell.type.labels = cell_type_labels, cell.state.labels = cell_state_lables_SCISSORS)

# In addition to the subclustering methods mentioned above, uses can also assign cell state labels based on their source of origin (using patient identifiers), when other subclustering methods are computationally intensive. This method is useful for cell types exhibiting clear intra-sample heterogeneity, such as malignant cells.

# for example, recluster Neurons based on their sampleID
cell_state_lables_manual <- cell_type_labels
Neurons_id <- which(pData(sc.eset)$cellType == 'Neurons')
cell_state_lables_manual[Neurons_id] <- paste0('Neurons_',pData(sc.eset)$SubjectName[Neurons_id])
table(cell_state_lables_manual)

# 3. Reference evaluation pipeline ----------------------------------------
# We recognize that cell state specification is an unsupervised problem and there is no definitive answer for the hyperparameter choices. Therefore, we have provided a reference evaluation pipeline to facilitate deconvolution performance evaluation and guide hyperparameter choices during reference construction.


# Import a BayesPrism object for deconvolution ----------------------------
# For users familiar with the BayesPrism working framework, InstaPrism offers a convenient way that directly accepts a Prism object derived from BayesPrism::new.prism() as input for fast deconvolution. Consider the following example from BayesPrism tutorial. The extdata folder associated can be downloaded from zenodo repository.

# 1. Built a Prism object following the BayesPrism tutorial ---------------
library(BayesPrism)
# the example data is downloaded from: https://github.com/Danko-Lab/BayesPrism/tree/main/tutorial.dat
load('../data/tutorial.gbm.rdata') 

sc.stat <- plot.scRNA.outlier(
  input = sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels = cell.type.labels,
  species = "hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw = TRUE #return the data used for plotting. 
  #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
)

bk.stat <- plot.bulk.outlier(
  bulk.input = bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID 
  sc.input = sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels = cell.type.labels,
  species = "hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw = TRUE
  #pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
)

# Filter outlier genes from scRNA-seq data
sc.dat.filtered <- cleanup.genes (input = sc.dat,
                                  input.type = "count.matrix",
                                  species = "hs", 
                                  gene.group = c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells = 5)
# Subset protein coding genes
sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
                                         gene.type = "protein_coding")
# construct a prism object
myPrism <- new.prism(
  reference=sc.dat.filtered.pc, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.state.labels,
  key="tumor",
  outlier.cut=0.01,
  outlier.fraction=0.1
)
# run BayesPrism
start.time <- Sys.time()
bp.res <- run.prism(prism = myPrism,
                    n.cores=20,
                    update.gibbs=T) # set update.gibbs = T, the default option
end.time <- Sys.time()
bp_running_time <- difftime(end.time, start.time,units = 'mins') %>% as.numeric()
save(bp.res,bp_running_time,file = '../out/object/bp.res.RData')

# Please note that running the above code took over 6 hours on our system (using 16 cores). Users may bypass the BayesPrism running process by directly loading our pre-processed results.
load('../out/object/bp.res.RData')

# With the default setting of BayesPrism (update.gibbs = TRUE), the BayesPrism::run.prism() function will performs two rounds of deconvolution: the first using an scRNA-based reference and the second with an updated reference. The estimated cell-type fractions are stored in bp.res@posterior.initial.cellType@theta and bp.res@posterior.theta_f@theta respectively.

# 2. Import the Prism object to InstaPrism --------------------------------
# The InstaPrism() function can directly accept a Prism object for deconvolution.
InstaPrism.res.initial <- InstaPrism(input_type = 'prism',
                                     prismObj = myPrism,
                                     n.core = 20)

# Note that InstaPrism() function only returns the posterior estimates using scRNA-based reference. The results will correspond to bp.res@posterior.initial.cellType returned by BayesPrism:
# compare cell type fraction estimates between two methods
ct.corr <- cor(t(InstaPrism.res.initial@Post.ini.ct@theta),bp.res@posterior.initial.cellType@theta) %>% round(1)

ComplexHeatmap::Heatmap(ct.corr,
                        show_row_dend = F,
                        show_column_dend = F,
                        column_title = 'BayesPrism',
                        row_title = 'InstaPrism',
                        name='correlation',
                        cell_fun = function(j, i, x, y, w, h, col) {grid::grid.text(ct.corr[i, j], x, y)})

# The cell-type specific expression are also nearly identical:
# tumor specific expression from InstaPrism
Z_tumor_ct_initial_InstaPrism = reconstruct_Z_ct_initial(InstaPrism_obj = InstaPrism.res.initial,cell.type.of.interest = 'tumor')

# tumor specific expression from BayesPrism
Z_tumor_ct_initial_bp = bp.res@posterior.initial.cellType@Z[,,'tumor']

# compare tumor specific expression from two methods
all.equal(t(Z_tumor_ct_initial_InstaPrism),Z_tumor_ct_initial_bp)
#> [1] "Mean relative difference: 0.003338094"

# 3. Get updated posterior with InstaPrism --------------------------------
# InstaPrism incorporates the capability of reference updates with the InstaPrism_update() function, matching the bp.res@posterior.theta_f returned by BayesPrism.
updated_obj = InstaPrism_update(InstaPrism.res.initial, 
                                bulk_Expr = t(bp.res@prism@mixture), 
                                key = 'tumor',
                                n.core = 16)
save(InstaPrism.res.initial,updated_obj,file = '../out/object/InstaPrism_res.RData')

# The updated posterior estimates of cell type fractions are still nearly identical betweent two methods
ct.corr.updated <- cor(t(updated_obj@theta),bp.res@posterior.theta_f@theta)
diag(ct.corr.updated)

# InstaPrism also supports deconvolution of gene expression using the updated reference, a feature not available with BayesPrism:
Z_updated <- get_Z_array(updated_obj) # the 3d expression array (sample by gene by cell type)
Z_tumor_ct_updated_InstaPrism <- reconstruct_Z_ct_updated(InstaPrism_updated_obj = updated_obj,cell.type.of.interest = 'tumor') # get the tumor specific expression only

# Additional notes on InstaPrism parameters -------------------------------
# So far we have demonstrated the primary usage of InstaPrism, mostly with the default parameters. In real practice, some of the default parameters can be tailored to meet the users’ specific requirements.

# n.iter: number of iterations in InstaPrism algorithm.
# n.iter is an important parameter for InstaPrism algorithm. InstaPrism updates the fraction estimates over iterations and n.iter determines how many iterations are performed. Under default setting, this value is set to max(100, number of cell.states ×
# 2). However, users can set it to larger values to improve convergence. Specifically, we included two additional parameters that help to verify the convergence status: verbose and convergence.plot.
# verbose: a logical variable to determine whether to display convergence status of the model.
# By setting verbose = T, InstaPrsim() will print the convergence status summary for both cell.states and cell.types.
# Specifically, the absolute difference in fraction estimates between the last two iterations is utilized as an indicator of convergence (abs_diff), with smaller values indicating convergence and higher confidence of the deconvolution results. Using example3.RData as a demonstration, by setting verbose = T, we get a summarized table showing the minimum, median and maximum abs_diff for all cell.states/types across the samples. Usually we consider abs_diff < 0.01 as convergence.

load('extdata/tutorial_example3/example3.RData')
ref <- refPrepare(scExpr_train,
                  cell.type.labels,
                  cell.state.labels)
InstaPrism.res = InstaPrism(bulk_Expr = simulated_bulk,refPhi_cs = ref,
                            n.iter = 100,n.core = 8,verbose = T)

# convergence.plot: a logical variable determining whether to visualize convergence status for cell types.
# Again, using example3.RData as a demonstration, first let’s consider the scenario whenn.iter is set too low. According to the convergence plot, fraction estimates for malignant cells did not converge. This means that we should set higher value of n.iter.
InstaPrism.res <- InstaPrism(bulk_Expr = simulated_bulk,
                             refPhi_cs = ref,
                             n.iter = 20,
                             n.core = 8,
                             verbose = T,
                             convergence.plot=T,
                             max_n_per_plot = 100)

# Now, let’s try using a higher value of n.iter. This time we get fraction estimation for all the cell types converge (with all the abs_diff less than 0.01).
InstaPrism.res <- InstaPrism(bulk_Expr = simulated_bulk,
                             refPhi_cs = ref,
                             n.iter = 100,
                             n.core = 8,
                             verbose = T,
                             convergence.plot=T,
                             max_n_per_plot = 100)

# n.core: number of cores to use for parallel programming. We recommend setting this parameter to higher value when multiple cores are available. InstaPrism functions that support parallel programming include: InstaPrism(), InstaPrism_update(), reconstruct_Z_ct_initial(), reconstruct_Z_ct_updated() and get_Z_array().

# Previous tutorial -------------------------------------------------------
# Our previous tutorial (associated with v0.1.5), which provides a more detailed comparison between InstaPrism and BayesPrism, can be accessed here.

