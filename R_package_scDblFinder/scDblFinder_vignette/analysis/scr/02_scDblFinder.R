# scDblFinder -------------------------------------------------------------
# Pierre-Luc Germain1
# 
# 1University and ETH Zürich
# 
# 15 April 2025
# Abstract
# An introduction to the scDblFinder method for fast and comprehensive doublet identification in single-cell data.
# 
# Package
# scDblFinder 1.22.0

# 1scDblFinder ------------------------------------------------------------
# The scDblFinder method combines the strengths of various doublet detection approaches, training an iterative classifier on the neighborhood of real cells and artificial doublets.
# scDblFinder() has two main modes of operation: cluster-based or not. Both perform quite well (see Germain et al., 2021). In general, we recommend the cluster-based approach in datasets with a very clear cluster structure, and the random approach in more complex datasets.

# 1.1Installation
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("scDblFinder")

# or, to get that latest developments:
# BiocManager::install("plger/scDblFinder")

# 1.2Usage
# The input of scDblFinder is an object sce of class SingleCellExperiment (empty drops having already been removed) containing at least the counts (assay ‘counts’). Alternatively, a simple count matrix can also be provided.
# Given an SCE object, scDblFinder (using the random approach) can be launched as follows :
set.seed(123)
# we create a dummy dataset; since it's small we set a higher doublet rate
sce <- mockDoubletSCE(dbl.rate=0.1, ngenes=300)
sce
sce@colData
# we run scDblFinder (providing the unusually high doublet rate)
sce <- scDblFinder(sce, dbr=0.1)
sce
sce@colData

# For 10x data, it is usually safe to leave the dbr empty, and it will be automatically estimated. (If using a chip other than the standard 10X, you might have to adjust it or the related dbr.per1k argument.
# scDblFinder will add a number of columns to the colData of sce prefixed with ‘scDblFinder’, the most important of which are:
# sce$scDblFinder.score : the final doublet score
# sce$scDblFinder.class : the classification (doublet or singlet)
# We can compare the calls with the truth in this toy example:
                                                                                                    
table(truth=sce$type, call=sce$scDblFinder.class)

# To use the cluster-based approach, one simply needs to additionally provide the clusters argument:
sce@colData$cluster %>% table()
sce <- scDblFinder(sce, clusters="cluster")
table(truth=sce$type, call=sce$scDblFinder.class)

# The clusters argument can be either a vector of cluster labels for each column of sce, a colData column of sce containing such labels, or TRUE. If clusters=TRUE, the fast clustering approach (see ?fastcluster) will be employed. If normalized expression (assay ‘logcounts’) and/or PCA (reducedDim ‘PCA’) are already present in the object, these will be used for the clustering step.
                                                                                                  
# 1.2.1Multiple samples
# If you have multiple samples (understood as different cell captures), then it is preferable to look for doublets separately for each sample (for multiplexed samples with cell hashes, this means for each batch). You can do this by simply providing a vector of the sample ids to the samples parameter of scDblFinder or, if these are stored in a column of colData, the name of the column. In this case, you might also consider multithreading it using the BPPARAM parameter (assuming you’ve got enough RAM!). For example:
library(BiocParallel)

# create some mock datasets
set.seed(123)
# we create a dummy dataset; since it's small we set a higher doublet rate
sce2 <- mockDoubletSCE(dbl.rate=0.1, ngenes=300)
sce3 <- mockDoubletSCE(dbl.rate=0.1, ngenes=300)

# add the sample id
sce2$sample_id <- "sce2"
sce3$sample_id <- "sce3"

# combine the two datasets
combined_sce <- cbind(sce2, sce3)
combined_sce
combined_sce@colData

combined_sce <- scDblFinder(combined_sce, samples="sample_id", BPPARAM=MulticoreParam(3))
combined_sce@colData
table(combined_sce$scDblFinder.class)

table(truth=combined_sce$type, call=combined_sce$scDblFinder.class)

# Note that if you are running multiple samples using the cluster-based approach (see below), clustering will be performed sample-wise. While this is typically not an issue for doublet identification, it means that the cluster labels (and putative origins of doublets) won’t match between samples. If you are interested in these, it is preferable to first cluster (for example using sce$cluster <- fastcluster(sce)) and then provide the clusters to scDblFinder, which will ensure concordant labels across samples.
# Of note, if you have very large differences in number of cells between samples the scores will not be directly comparable. We are working on improving this, but in the meantime it would be preferable to stratify similar samples and threshold the sets separately.

# 1.3Description of the method --------------------------------------------
# Wrapped in the scDblFinder function are the following steps:
# 1.3.1Splitting captures
# Doublets can only arise within a given sample or capture, and for this reason are better sought independently for each sample, which also speeds up the analysis. If the samples argument is given, scDblFinder will use it to split the cells into samples/captures, and process each of them in parallel if the BPPARAM argument is given. Depending on the multiSampleMode argument, the classifier can be trained globally, with thresholds optimized on a per-sample basis; however we did not see an improvement in doing so, and therefore by default each sample is treated separately to maximize robustness to technical differences.
# 
# If your samples are multiplexed, i.e. the different samples are mixed in different batches, then the batches should be what you provide to this argument.
# 
# 1.3.2Reducing and clustering the data
# The analysis can be considerably sped up, at little if any cost in accuracy, by reducing the dataset to only the top expressed genes (controlled by the nfeatures argument).
# 
# Then, depending on the clusters argument, an eventual PCA and clustering (using the internal fastcluster function) will be performed. The rationale for the cluster-based approach is that homotypic doublets are nearly impossible to distinguish on the basis of their transcriptome, and therefore that creating that kind of doublets is a waste of computational resources that can moreover mislead the classifier into flagging singlets. An alternative approach, however, is to generate doublets randomly (setting clusters to FALSE or NULL), and use the iterative approach (see below) to exclude also unidentifiable artificial doublets from the training.
# 
# 1.3.3Generating artificial doublets
# Depending on the clusters and propRandom arguments, artificial doublets will be generated by combining random cells and/or pairs of non-identical clusters (this can be performed manually using the getArtificialDoublets function). A proportion of the doublets will simply use the sum of counts of the composing cells, while the rest will undergo a library size adjustment and poisson resampling.
# 
# 1.3.4Examining the k-nearest neighbors (kNN) of each cell
# A new PCA is performed on the combination of real and artificial cells, from which a kNN network is generated. Using this kNN, a number of parameters are gathered for each cell, such as the proportion of doublets (i.e. artificial doublets or known doublets provided through the knownDoublets argument, if given) among the KNN, ratio of the distances to the nearest doublet and nearest non-doublet, etc. Several of this features are reported in the output with the ‘scDblFinder.’ prefix, e.g.:
#   
#   distanceToNearest : distance to the nearest cell (real or artificial)
# ratio : the proportion of the KNN that are doublets. (If more than one value of k is given, the various ratios will be used during classification and will be reported)
# weighted : the proportion of the KNN that are doublets, weighted by their distance (useful for isolated cells)
# 1.3.5Training a classifier
# Unless the score argument is set to ‘weighted’ or ‘ratio’ (in which case the aforementioned ratio is directly used as a doublet score), scDblFinder then uses gradient boosted trees trained on the kNN-derived properties along with a few additional features (e.g. library size, number of non-zero features, and an estimate of the difficultly of detecting artificial doublets in the cell’s neighborhood, a variant of the cxds score from the scds, etc.) to distinguish doublets (either artificial or given) from other cells, and assigns a score on this basis.
# 
# One problem of using a classifier for this task is that some of the real cells (the actual doublets) are mislabeled as singlet, so to speak. scDblFinder therefore iteratively retrains the classifier, each time excluding from the training the (real) cells called as doublets in the previous step (as well as unidentifiable artificial doublets). The number of steps being controlled by the iter parameter (in our experience, 2 or 3 is optimal).
# 
# This score is available in the output, in the scDblFinder.score colData column, and can be interpreted as a probability. If the data is multi-sample, a single model is trained for all samples.
# 
# 1.3.6Thresholding
# Rather than thresholding on some arbitrary cutoff of the score, scDblFinder uses the expected number of doublets in combination to the misclassification rate to establish a threshold. Unless it is manually given through the dbr argument, the expected doublet rate is first estimated (see below). If samples were specified, and if the dbr is automatically calculated, thresholding is performed separately across samples.
# 
# Thresholding then tries to simultaneously minimize: 1) the classification error (in terms of the proportion of known doublets below the threshold) and 2) the deviation from the expected number of doublets among real cells (as a ratio of the total number of expected doublets within the range determined by dbr.sd, and adjusted for homotypic doublets). This means that, if you have no idea about the doublet rate, setting dbr.sd=1 will make the threshold depend entirely on the misclassification rate.
# 
# 1.3.7Doublet origins and enrichments
# If artificial doublets are generated between clusters, it is sometimes possible to call the most likely origin (in terms of the combination of clusters) of a given putative real doublet. We observed that at least one of the two composing cell is typically recognized, but that both are seldom correctly recognized, owing to the sometimes small relative contribution of one of the two original cells. This information is provided through the scDblFinder.mostLikelyOrigin column of the output (and the scDblFinder.originAmbiguous column indicates whether this origin is ambiguous or rather clear). This, in turn, allows us to identify enrichment over expectation for specific kinds of doublets. Some statistics on each combination of clusters are saved in metadata(sce)$scDblFinder.stats, and the plotDoubletMap function can be used to visualize enrichments. In addition, two frameworks are offered for testing the significance of enrichments:
# sample usage of the function:

library(SingleCellExperiment)
library(scater)
library(scran) 
library(circlize)
library(ComplexHeatmap)

set.seed(123) # for reproducibility
sce <- mockDoubletSCE(dbl.rate = 0.1, ngenes = 500)

# b. Normalization
sce <- logNormCounts(sce)

# c. Feature selection (optional but often good for UMAP)
# Identify highly variable genes
hvgs <- getTopHVGs(sce, n = 2000) # Select top 2000 HVGs

# d. Dimensionality Reduction (e.g., PCA)
sce <- runPCA(sce, subset_row = hvgs)

# e. Clustering (optional, but good for context on UMAP)
# Use a graph-based clustering method
snn_graph <- buildSNNGraph(sce, k = 10, use.dimred = "PCA")
sce$cluster <- factor(igraph::cluster_louvain(snn_graph)$membership)

# f. Compute UMAP
# Ensure you specify the reducedDim "PCA" that you just computed
sce <- runUMAP(sce, dimred = "PCA")

# 3. Run scDblFinder
# This will add 'scDblFinder.class' and 'scDblFinder.score' to colData(sce)
sce <- scDblFinder(sce,clusters = "cluster")

# 4. Sample Usage of plotDoubletMap
plotDoubletMap(sce)

# The clusterStickiness function tests whether each cluster forms more doublet than would be expected given its abundance, by default using a single quasi-binomial model fitted across all doublet types.
# The doubletPairwiseEnrichment function separately tests whether each specific doublet type (i.e. combination of clusters) is more abundant than expected, by default using a poisson model.

# 1.4Some important parameters --------------------------------------------
# scDblFinder has a fair number of parameters governing the preprocessing, generation of doublets, classification, etc. (see ?scDblFinder). Here we describe just a few of the most important ones.
# 
# 1.4.1Expected proportion of doublets
# The expected proportion of doublets has no impact on the density of artificial doublets in the neighborhood, but impacts the classifier’s score and, especially, where the cutoff will be placed. It is specified through the dbr parameter, as well as dbr.per1k (which specifies, if dbr is omitted, the rate per thousands cells from which to estimate it). In addition, the dbr.sd parameter specifies a +/- range around dbr within which the deviation from dbr will be considered null.
# 
# For most platforms, the more cells you capture the higher the chance of creating a doublet. For standard 10X data, the 10X documentation indicates a doublet rate of roughly 0.8% per 1000 cells captured, which is the default value of dbr.per1k. This means that unless dbr is manually set, with 5000 cells, (0.008*5)*5000 = 200 doublets are expected, and the default expected doublet rate will be set to this value (with a default standard deviation of 0.015). Note however that different protocols may vary in the expected proportion of doublets. For example, the high-throughput (HT) 10X kit has an expected doublet rate of half the standard, i.e. 0.4% per 1000 cells, so if using that kit, set dbr.per1k=0.004.
# 
# Also note that strictly speaking, the proportion of doublets depends more on the number of cells inputted than that recovered. If your recovery rate was lower than expected, you might observe a higher doublet rate (see the too-many doublets section below).
# 
# The impact of the expected doublet rate on the thresholding will depend on how hard the classification task is: if it is easy, the called doublets will not depend much on the expected rate. If you are unsure about the doublet rate, you might consider increasing dbr.sd: with a high value (e.g. 1), the thresholding will be entirely based on the misclassification error (without any assumption about an expected doublet rate).
# 
# 1.4.2Number of artificial doublets
# The number of artificial doublets can be set through the artificialDoublets parameter. Using more artificial doublets leads to a better sampling of the possible mixtures of cells, but increases memory and runtime, and can skew the scores, in extreme cases leading to difficulties in setting a threshold for being called as a doublet (see this issue for a discussion).
# 
# By default, scDblFinder will generate roughly as many artificial doublets as there are cells, which is usually appropriate. However, for very small datasets this could represent an undersampling of the mixing space and hence lead to lower detection accuracy. For this reason, a hard minimum number of artificial doublet is set. This will tend to improve accuracy for small datasets, but the scores will be skewed towards 1, possibly making a separation difficult. If you are in such a situation and your histogram of scores does not show a bimodality, consider manually setting the artificialDoublets parameter to something closer to your actual number of cells.
# 
# 
# 
# 
# 1.5Frequently-asked questions
# 1.5.1I’m getting way too many doublets called - what’s going on?
#   Then you most likely have a wrong doublet rate. If you did not provide it (dbr argument), the doublet rate will be calculated automatically using expected doublet rates from 10x, meaning that the more cells captured, the higher the doublet rates. If you have reasons to think that this is not applicable to your data, set the dbr manually.
# 
# The most common cause for an unexpectedly large proportion of doublets is if you have a multi-sample dataset and did not split by samples. scDblFinder will think that the data is a single capture with loads of cells, and hence with a very high doublet rate. Splitting by sample should solve the issue.
# 
# Also note that, although 10X-like data tends to have roughly 1% per 1000 cells captured, the determining factor for doublet formation is the number of cells inserted into the machine. If for some reason your recovery rate is lower than expected, you might have a higher doublet rate than you’d expect from the captured and called cells (in other words, it would be preferable to say that the doublet rate is roughly 0.6% per 1000 cells put into the machine, where 0.6 is the recovery rate). In such circumstances, scDblFinder typically sets the thresholds correctly nevertheless. This is because the thresholding tries to minimize both the deviation from the expected number of doublets and the misclassification (i.e. of artificial doublets), meaning that the effective (i.e. final) doublet rate will differ from the given one. scDblFinder also considers false positives to be less problematic than false negatives. You can reduce to some degree the deviation from the input doublet rate by setting dbr.sd=0.
# 
# Finally, note that version (1.20.0) initially shipped with the current Bioconductor release (3.20) version included a wrong default doublet rate (dbr.per1k) argument (it was 0.08 instead of 0.008). This was subsequently fixed in version 1.20.2, but if you installed before that you might need to update the package.
# 
# 1.5.2Should I use the cluster-based doublet generation or not?
#   Both approaches perform very similarly overall in benchmarks (see the scDblFinder paper). If your data is very clearly segregated into clusters, or if you are interested in the origin of the doublets, the cluster-based approach is preferable. This will also enable a more accurate accounting of homotypic doublets, and therefore a slightly better thresholding. Otherwise, and especially if your data does not segregate very clearly into clusters, the random approach (e.g. clusters=FALSE, the default) is preferable.
# 
# 1.5.3The clusters don’t make any sense!
#   If you ran scDblFinder on a multi-sample dataset and did not provide the cluster labels, then the labels are sample-specific (meaning that label ‘1’ in one sample might have nothing to do with label ‘1’ in another), and plotting them on a tSNE will look like they do not make sense. For this reason, when running multiple samples we recommend to first cluster all samples together (for example using sce$cluster <- fastcluster(sce)) and then provide the clusters to scDblFinder.
# 
# 1.5.4‘Size factors should be positive’ error
# You will get this error if you have some cells that have zero reads (or a very low read count, leading to zero after feature selection). After filtering out these cells the error should go away.
# 
# 1.5.5Identifying homotypic doublets
# Like other similar tools, scDblFinder focuses on identifying heterotypic doublets (formed by different cell types), and has only a low performance in identifying homotypic doublets (see this preprint). This can lead to disagreements with doublets called using cell hashes or SNPs in multiplexed samples, which capture both types of doublets similarly (and can miss intra-sample heterotypic doublets, especially if the multiplexing is low). This is why we treat these approaches as complementary.
# 
# However, should you for some reason try to identify also homotypic doublets with scDblFinder, be sure to not to use the cluster-based approach, and to set removeUnidentifiable=FALSE. Otherwise, scDblFinder removes artificial doublets likely to be homotypic from training, therefore focusing the task on heterotypic doublets, but at the expense ot homotypic ones (which are typically deemed relatively harmless).
# 
# 1.5.6What is a sample exactly? Usage with barcoded and 10X Flex data.
# As indicated above, the samples argument should be used to indicate different captures. For multiplexed samples, this is expected to be the batch of cells processed together, rather than the actual samples.
# 
# In highly multiplexed datasets such as produced by the 10X Flex kit (especially 16-plex), this can cause two kinds of problems. First, the whole logic of the Flex approach is that inter-sample doublets can be resolved into separate cells, and while a large number of unresolvable intra-sample doublets will remain (see Howitt et al., 2024), the expected remaining doublet rate will not be the same as for classical 10X experiment. For this reason, we recommend to set a higher dbr.sd in such circumstances, e.g. dbr.sd=1 to base the thresholding entirely on the classification accuracy.
# 
# Another, more practical problem is that, with such kits, the very large number of cells in a single capture might translante into very large computational demands when running scDblFinder. To circumvent such problem, one can split a batch of cells into more decently-sized chunks and process the chunks separately, so long as each chunk is representative of the whole batch in terms of cell heterogeneity.
# 
# 1.5.7How can I make this reproducible?
#   Because it relies on the partly random generation of artificial doublets, running scDblFinder multiple times on the same data will yield slightly different results. You can ensure reproducibility using set.seed(), however this will not be sufficient when processing multiple samples (i.e. using the samples argument – even without multithreading!). In such case, the seed needs to be passed to the BPPARAMs:

set.seed(123)
# we create a dummy dataset; since it's small we set a higher doublet rate
sce2 <- mockDoubletSCE(dbl.rate=0.1, ngenes=300)
sce3 <- mockDoubletSCE(dbl.rate=0.1, ngenes=300)

# add the sample id
sce2$sample_id <- "sce2"
sce3$sample_id <- "sce3"

# combine the two datasets
combined_sce <- cbind(sce2, sce3)

# run whithout seed
combined_sce0 <- scDblFinder(combined_sce, clusters="cluster", samples="sample_id", BPPARAM=MulticoreParam(10))
combined_sce02 <- scDblFinder(combined_sce, clusters="cluster", samples="sample_id", BPPARAM=MulticoreParam(10))

# confirm the result are different
all.equal(combined_sce0$scDblFinder.score,combined_sce02$scDblFinder.score)

# run using the seed
bp <- MulticoreParam(10, RNGseed=1234)
combined_sce <- scDblFinder(combined_sce, clusters="cluster", samples="sample_id", BPPARAM=bp)
combined_sce2 <- scDblFinder(combined_sce, clusters="cluster", samples="sample_id", BPPARAM=bp)

# confirm the result is the same
all.equal(combined_sce$scDblFinder.score,combined_sce2$scDblFinder.score)


# Similarly, when processing the samples serially, use SerialParam(RNGseed = seed).
# 
# (Note that in BiocParallel versions <1.28, one had in addition to explicitly start the cluster before the run using bpstart(bp), and then bpstop(bp) after scDblFinder.)
# 
# As a final note: when running scDblFinder twice on the same data with different random seeds, the scores will be highly correlated, but some cells will be called as doublets (with a high score) in only one of the runs (e.g. see this issue). There are good reasons to believe that these are homotypic doublets (if doublets at all), and if you worry chiefly about hetertypic doublets, you may concentrate on those that are reprocibly called across runs.
# 
# 1.5.8Can I use this in combination with Seurat or other tools?
#   If the input SCE already contains a logcounts assay or a reducedDim slot named ‘PCA’, scDblFinder will use them for the clustering step. In addition, a clustering can be manually given using the clusters argument of scDblFinder(). In this way, seurat clustering could for instance be used to create the artificial doublets (see ?Seurat::as.SingleCellExperiment.Seurat for conversion to SCE). For example, assuming as Seurat object se, the following could be done:

# simulate the creation of a seurat object
set.seed(123) # for reproducibility
sce <- mockDoubletSCE(dbl.rate = 0.1, ngenes = 500)

# b. Normalization
sce <- logNormCounts(sce)

# c. Feature selection (optional but often good for UMAP)
# Identify highly variable genes
hvgs <- getTopHVGs(sce, n = 2000) # Select top 2000 HVGs

# d. Dimensionality Reduction (e.g., PCA)
sce <- runPCA(sce, subset_row = hvgs)

# e. Clustering (optional, but good for context on UMAP)
# Use a graph-based clustering method
snn_graph <- buildSNNGraph(sce, k = 10, use.dimred = "PCA")
sce$cluster <- factor(igraph::cluster_louvain(snn_graph)$membership)

# f. Compute UMAP
# Ensure you specify the reducedDim "PCA" that you just computed
sce <- runUMAP(sce, dimred = "PCA")

# convert to seurat
se <- Seurat::as.Seurat(sce)
# define the cluster identity
Idents(se) <- se$cluster

se

# pull the raw counst and the cluster identiry of the object
set.seed(123)
sce_se <- scDblFinder(GetAssayData(se, slot="counts"), clusters=Idents(se))

# port the resulting scores back to the Seurat object:
se$scDblFinder.score <- sce_se$scDblFinder.score

# is this different from runnign scDblFinder on the original object?
set.seed(123)
sce_test <- scDblFinder(sce, clusters="cluster")

# compare the results of the scores
all.equal(sce_test$scDblFinder.score,unname(se$scDblFinder.score))


# After artificial doublets generation, the counts of real and artificial cells must then be reprocessed (i.e. normalization and PCA) together, which is performed internally using scater. If you wish this step to be performed differently, you may provide your own function for doing so (see the processing argument in ?scDblFinder). We note, however, that the impact of variations of this step on doublet detection is rather mild. In fact, not performing any normalization at all for instance decreases doublet identification accuracy, but by rather little.
# 
# For example, the following code would enable the internal use of sctransform:
  
# assuming `x` is the count matrix:
nfeatures <- 1000
sce <- SingleCellExperiment(list(counts=x))
# sctransform on real cells:
vst1 <- sctransform::vst(counts(sce), n_cells=min(ncol(sce),5000), verbosity=0)
sce <- sce[row.names(vst1$y),]
logcounts(sce) <- vst1$y
hvg <- row.names(sce)[head(order(vst1$gene_attr$residual_variance, decreasing=TRUE), nfeatures)]

# define a processing function that scDblFinder will use on the real+artificial doublets;
# the input should be a count matrix and the number of dimensions, and the output a PCA matrix

myfun <- function(e, dims){
  # we use the thetas calculated from the first vst on real cells
  e <- e[intersect(row.names(e), row.names(vst1$model_pars_fit[which(!is.na(vst1$model_pars_fit[,"theta"])),])),]
  vst2 <- sctransform::vst(e, n_cells=min(ncol(e),5000), method="nb_theta_given", 
                           theta_given=vst1$model_pars_fit[row.names(e),"theta"],
                           min_cells=1L, verbosity=0)
  scater::calculatePCA(vst2$y, ncomponents=dims)
}

sce <- scDblFinder(sce, processing=myfun, nfeatures=hvg)
# Note however that this did not generally lead to improved performance – but rather decreased on most benchmark datasets, in fact (see comparison in this issue).

# 1.5.9How can I call scDblFinder from the command line?
#   Here would be an example of how to call scDblFinder (in cluster mode) from the command line and save the results to a csv:
#   
#   Rscript -e '
# library(scDblFinder)
# set.seed(123) # for reproducibility
# e <- Matrix::readMM("matrix.mtx.gz")
# colnames(e) <- readLines("barcodes.tsv.gz")
# res <- scDblFinder(e, cluster=TRUE)
# res <- cbind(barcode=colnames(res),
#              colData(res)[,grep("scDblFinder",colnames(colData(res)))])
# write.table(res, "output.csv", row.names=FALSE, quote=FALSE)
# '
# 1.5.10Can this be used with scATACseq data?
#   Yes, see the scATAC vignette specifically on this topic.
# 
# 1.5.11Should I run QC cell filtering before or after doublet detection?
#   The input to scDblFinder should not include empty droplets, and it might be necessary to remove cells with a very low coverage (e.g. <200 or 500 reads) to avoid errors. Further quality filtering should be performed downstream of doublet detection, for two reasons: 1. the default expected doublet rate is calculated on the basis of the cells given, and if you excluded a lot of cells as low quality, scDblFinder might think that the doublet rate should be lower than it is. 2. kicking out all low quality cells first might hamper our ability to detect doublets that are formed by the combination of a good quality cell with a low-quality one. This being said, these are mostly theoretical grounds, and unless your QC filtering is very stringent (and it shouldn’t be!), it’s unlikely to make a big difference.
# 
# 1.5.11.1What about ambiant RNA decontamination?
#   Contamination by ambiant RNA has emerged as an important confounder in single-cell (and especially single-nuclei) RNAseq data, which prompts the question of whether that should be run prior or after doublet detection. Unfortunately, we do not currently have good evidence pointing in either direction, and arguments can be made for both. Low-quality doublets, or doublets from an experiment with a large dominant celltype, can easily look like contamination, and likewise a high amount of contamination can easily look like a doublet because it includes RNA from other cell types. There is a possibility that a decontamination package sees an actual doublet as contamination, and attempts to clean it, which it will necessarily do imperfectly (because while the decontamination is a mixture of all cells, a doublet isn’t), but perhaps sufficiently so that it can’t be accurately detected as a doublet anymore. This would therefore be an argument for running doublet calling first. However, it’s also possible that decontamination, because it makes the cells cleaner, makes the doublet detection task easier.
# 
# 1.5.12Can I combine this method with others?
#   Of course it is always possible to run multiple methods and combine the results. In our benchmark, the combination of scDblFinder with DoubletFinder, for instance, did yield an improvement in most (though not all) datasets (see the results here), although of a small magnitude. The simplest way is to do an average of the scores (assuming that the scores are on a similar scale, and that a higher score has the same interpretation across methods), which for instance gave similar results to using a Fisher p-value combination on 1-score (interpreted as a probability).