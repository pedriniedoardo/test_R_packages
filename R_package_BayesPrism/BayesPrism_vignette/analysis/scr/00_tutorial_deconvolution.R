# AIM ---------------------------------------------------------------------
# run the tutorial analysis for the BayesPrism package
# usethis::create_github_token()
# gitcreds::gitcreds_set()

# libraries ---------------------------------------------------------------
library(tidyverse)
library(NMF)
library(BayesPrism)

# Load the dataset --------------------------------------------------------
# The rdata file is stored at "your_git_repository/tutorial.dat/tutorial.gbm.rdata"

load("../../ref/BayesPrism-main/tutorial.dat/tutorial.gbm.rdata")
ls()

# The rdata file contains four objects required for running BayesPrism.
str(bk.dat)
str(sc.dat)

str(cell.state.labels)
str(cell.type.labels)

# bk.dat: The sample-by-gene raw count matrix of bulk RNA-seq expression.
# rownames are bulk sample IDs, while colnames are gene names/IDs.
dim(bk.dat)
head(rownames(bk.dat))
head(colnames(bk.dat))

# sc.dat: The cell-by-gene raw count matrix of single cell RNA-seq expression.
# rownames are cell IDs, while colnames are gene names/IDs. sc.dat should be a dense matrix.
# If your sc.dat is a sparse matrix (dgCMatrix), you should convert it to a dense matrix, i.e., sc.dat <- as.matrix(sc.dat)
dim(sc.dat)
head(rownames(sc.dat))
head(colnames(sc.dat))

# Note that BayesPrism will perform deconvolution over genes shared between the mixture and scRNA-seq reference, i.e., by intersecting colnames(mixture) and colnames(reference). Therefore, please make sure to use consistent gene annotations (TCGA RNA-seq uses GENCODE v22).

# We recommend the use of unnormalized and untransformed count data. When raw count is not available, linear normalization, such as TPM, RPM, RPKM, FPKM, is also acceptable, as BayesPrism is robust to linear multiplicative difference between the reference and mixture. Ideally, if using normalized data, it is the best to supply reference and bulk transformed by the same method. log transformation should be avoided.

# cell.type.labels is a character vector of the same length as nrow(sc.dat) to denote the cell type of each cell in the reference.
sort(table(cell.type.labels))

# cell.state.labels is a character vector of the same length as nrow(sc.dat) to denote the cell state of each cell in the reference. In our example, cell states of malignant cells were obtained by sub-clustering the malignant cells from each patient, and cell states of myeloid cells were obtained by clustering myeloid cells from all patients. We define multiple cell states for these two cell types, as they contain substantial heterogeneity while also having sufficient number of cells for sub-clustering.
sort(table(cell.state.labels))

# check the summary table for the two vectors
table(cbind.data.frame(cell.state.labels, cell.type.labels))

# Please make sure that all cell states contain a reasonable number of cells, e.g. >20 or >50, so that their profile can be represented accurately.

# What to supply for cell.state.labels and cell.type.labels? The definition of cell type and cell state can be somewhat arbitrary (similar to the issue of assigning cell types for scRNA-seq) and depends on the question of interest. Their definitions depend on the granularity we aim at and the confidence of the cell.type.labels in scRNA-seq data. Usually, a good rule of thumb is as follows.

# 1) Define cell types as the cluster of cells having a sufficient number of significantly differentially expressed genes than other cell types, e.g., greater than 50 or even 100. For clusters that are too similar in transcription, we recommend treating them as cell states, which will be summed up before the final Gibbs sampling. Therefore, cell states are often suitable for cells that form a continuum on the phenotypic manifold rather than distinct clusters.

# 2) Define multiple cell states for cell types of significant heterogeneity, such as malignant cells, and of interest to deconvolve their transcription.


# QC of cell type and state labels ----------------------------------------
# We recommend first plotting the pairwise correlation matrix between cell states and between cell types. This will give us a sense of their quality. In cases where cell types/states are not represented by sufficient amount of information (low cell count and/or low library size), the low-quality cell types/states tend to cluster together. Users may re-cluster the data at higher granularity, or merge those cell types/states with the most similar cell types/states, or remove them (if re-clustering and merging is not appropriate).

plot.cor.phi(input = sc.dat,
             input.labels = cell.state.labels,
             title="cell state correlation",
             #specify pdf.prefix if need to output to pdf
             #pdf.prefix="gbm.cor.cs", 
             cexRow=0.2, cexCol=0.2,
             margins=c(2,2))


plot.cor.phi(input = sc.dat,
             input.labels = cell.type.labels,
             title="cell type correlation",
             #specify pdf.prefix if need to output to pdf
             #pdf.prefix="gbm.cor.ct",
             cexRow=0.5, cexCol=0.5
             )

# Filter outlier genes ----------------------------------------------------
# Gene expressed at high magnitude, such as ribosomal protein genes and mitochondrial genes, may dominate the distribution and bias the inference. These genes are often not informative in distinguishing cell types and can be a source of large spurious variance. As a result, they can be detrimental to deconvolution. We recommend the removal of these genes.
# Users may visualize the distribution of outlier genes from scRNA-seq reference. We compute the mean expression and of each gene across all cell types, and their cell type specificity scores.
# Visualize and determine outlier genes from scRNA-seq data

sc.stat <- plot.scRNA.outlier(
  input = sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels = cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE #return the data used for plotting. 
  #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
)

# As shown by the plot, ribosomal protein genes often show high mean expression and low cell type specificity scores.
# Users may also subset genes from sc.dat based on the statistics outputted by the function if needed.

head(sc.stat)  

# sc.stat shows the log of normalized mean expression (x-axis) and the maximum specificity (y-axis) of each gene, and if each gene belongs to a potential outlier category. other_Rb are a group of curated genes mostly consists of ribosomal pseudo-genes.

# Similarly, we may also visualize outlier genes from bulk RNA-seq. We compute the mean expression and of each gene across all cell types. As we do not have cell type level information from bulk data, we compute cell type specificity score from scRNA-seq, same as above.
# Visualize outlier genes in bulk RNA-seq
bk.stat <- plot.bulk.outlier(
  bulk.input = bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID 
  sc.input = sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels = cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE
  #pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
)

# check statistics from bulk RNA-seq data
head(bk.stat)

# Filter outlier genes from scRNA-seq data
# Next, we remove the genes from selected groups. Note that when sex is not identical between the reference and mixture, we recommend excluding genes from chrX and chrY. We also remove lowly transcribed genes, as the measurement of transcription of these genes tend to be noise-prone. Removal of these genes can also speed up computation.

sc.dat.filtered <- cleanup.genes(input = sc.dat,
                                 input.type = "count.matrix",
                                 species = "hs",
                                 gene.group = c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                 exp.cells = 5)

dim(sc.dat.filtered)
dim(sc.dat)

# Next, we check the concordance of gene expression for different types of genes. As bulk and single cell data are usually collected by different experimental protocols, they may have different sensitivity to different types of genes.
# note this function only works for human data. For other species, you are advised to make plots by yourself.
plot.bulk.vs.sc(sc.input = sc.dat.filtered,
                bulk.input = bk.dat
                #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
                )

# We observe that protein coding genes are the most concordant group between two assays. To reduce batch effects and speed up computation, we perform deconvolution on protein coding genes. We have also tried to run BayesPrism on all genes. The results were similar.

# Subset protein coding genes.
sc.dat.filtered.pc <-  select.gene.type(sc.dat.filtered,
                                        gene.type = "protein_coding")
str(sc.dat.filtered.pc)

# Optionally, in cases where cell types are defined in a way that some of them show very similar transcription or severe batch effects exist, e.g., reference and mixture are from very different assays (ribo-depleted RNA-seq vs poly-A tail RNA-seq or PRO-seq (nascent RNA-seq) vs RNA-seq (steady state RNA)), selecting signature genes can be beneficial. This is because the selection of signature genes can enrich for genes informative for deconvolution while reducing the impact of noise caused by technical batch effects.
# We provide a function for selecting genes by performing differential expression using pair-wise t-test between cell states from different cell types. Other differential expression analysis can also be used.


# Select marker genes (Optional)
# performing pair-wise t test for cell states from different cell types

diff.exp.stat <- get.exp.stat(sc.dat = sc.dat[,colSums(sc.dat>0)>3],# filter genes to reduce memory use
                              cell.type.labels = cell.type.labels,
                              cell.state.labels = cell.state.labels,
                              pseudo.count = 0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                              cell.count.cutoff = 50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                              n.cores = 20 #number of threads
                              )

str(diff.exp.stat)

# Ideally, we would like to have sufficient number of genes selected for each cell type (>50). If not, users may lower the cutoff.

# To subset our count matrix over the signature genes, do

sc.dat.filtered.pc.sig <- select.marker(sc.dat = sc.dat.filtered.pc,
                                        stat = diff.exp.stat,
                                        pval.max = 0.01,
                                        lfc.min = 0.1)

dim(sc.dat.filtered.pc.sig)
sc.dat.filtered.pc.sig[1:10,1:10]

# To run BayesPrism using the signature genes, use reference=sc.dat.filtered.sig when call new.prism (see below).


# Construct a prism object. -----------------------------------------------
# A prism object contains all data required for running BayesPrism, namely, a scRNA-seq reference matrix, the cell type and cell state labels of each row of reference, and the mixture matrix for bulk RNA-seq. 

# When using scRNA-seq count matrix as the input (recommended), user needs to specify input.type = "count.matrix". The other option for input.type is "GEP" (gene expression profile) which is a cell state by gene matrix. This option is used when using reference derived from other assays, such as sorted bulk data.

# The parameter key is a character in cell.type.labels that corresponds to the malignant cell type. Set to NULL if there are no malignant cells or the malignant cells between reference and mixture are from matched samples, in which case all cell types will be treated equally.

myPrism <- new.prism(
  reference = sc.dat.filtered.pc, 
  mixture = bk.dat,
  input.type = "count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.state.labels,
  key="tumor",
  outlier.cut=0.01,
  outlier.fraction=0.1,
  )

# Note that outlier.cut and outlier.fraction=0.1 filter genes in X whose expression fraction is greater than outlier.cut (Default=0.01) in more than outlier.fraction (Default=0.1) of bulk data. 
# Typically for dataset with reasonable quality control, very few genes will be filtered. 
# Removal of outlier genes will ensure that the inference will not be dominated by outliers, which sometimes may be resulted from poor QC in mapping.

# Run BayesPrism. ---------------------------------------------------------
# Next, we start the run of BayesPrism.
# Parameters to control Gibbs sampling and optimization can be specified using gibbs.control and opt.control. Do ?run.prism for details. We recommend the use of default parameters.

bp.res <- run.prism(prism = myPrism, n.cores=50)

# Report the summary statistics.

bp.res

saveRDS(bp.res,"../out/object/bp.res.rds")

# Alternatively, if user would like to run step 4 separately from 1-3, one may do the following:
  
# not run
# bp.res.initial <- run.prism(prism = myPrism, n.cores=50, update.gibbs=FALSE)
# bp.res.update <- update.theta (bp = bp.res.initial)


# Extract results ---------------------------------------------------------
# Understand the results BayesPrism is expected to generate the following results: 
# The output from bp.res is an S4 object of the class "BayesPrism". Let’s take a look at what’s inside:
slotNames(bp.res)

# "prism" is the input prism.
# "posterior.initial.cellState" is the result of step 2.
# "posterior.initial.cellType" is the result of step 3.
# "reference.update" is the updated reference ψ
# .
# "posterior.theta_f" is the result of step 4.
# "control_param" contains the parameters to run BayesPrism.

# We provide utility functions to extract deconvolved cell type fraction and expression from the output of run.prism.

# extract posterior mean of cell type fraction theta
theta <- get.fraction(bp=bp.res,
                      which.theta="final",
                      state.or.type="type")

head(theta)

# extract coefficient of variation (CV) of cell type fraction
theta.cv <- bp.res@posterior.theta_f@theta.cv

head(theta.cv)

# theta.cv quantifies how concentrated the posterior distribution is. Higher theta value is associated with lower CV. When performing rank-ordered statistics, such as Spearman’s rank correlation, users may mask theta (e.g. clip at zero for small theta) where the CV is above some number, say > 0.1. Lower sequencing depth, which is the case for spatial transcriptomics, is generally associated with higher CV. A reasonable cutoff of CV for deconvolving Visium data may be ~0.5.

# extract posterior mean of cell type-specific gene expression count matrix Z  
Z.tumor <- get.exp(bp=bp.res,
                   state.or.type="type",
                   cell.name="tumor")

head(t(Z.tumor[1:5,]))

# try to extract the bulk expression for another annotaiton
Z.endo <- get.exp(bp=bp.res,
                  state.or.type="type",
                  cell.name="endothelial")

head(t(Z.endo[1:5,]))

# save the result
save(bp.res, file="bp.res.rdata")

# Downstream analysis -----------------------------------------------------
# Potential downstream analysis can be performed using theta and Z are:
# Clustering bulk samples by theta or Z (Z can be normalized by vst(round(t(Z.tumor))), using the vst function from the DESeq2 package.)
# Computing z-scores of signature genes for Z of the cell type of interest.
# Survival analysis and correlation with other clinical covariates using theta and Z.
# Correlating Z (after normalization using vst or from bp.res@reference.update@psi_mal) with theta to understand how gene expression of each gene (in malignant cells) correlates with the cell type fraction of non-malignant cells in tumor micro-environment, followed by gene set enrichment analysis (as done in BayesPrism paper).
# Embedding learning of malignant gene expression (see “Tutorial: embedding learning of malignant cell expression using BayesPrism”)

# example -----------------------------------------------------------------
# how to use the deconvoluted values for analysis
# Differential analysis on deconvolved cell type expression
# https://github.com/Danko-Lab/BayesPrism/issues/121
# Hello. To use DEseq2, you should 1) use raw count as input for BayesPrism, and 2) round up the Z matrix to the nearest integer and use it as input for DESeq2 (without any normalization step).

# Can BayesPrism be used to deconvolute normal and disease muscle tissues to different cell types?
# https://github.com/Danko-Lab/BayesPrism/issues/76
# Yes. BayesPrism works for both tumor and non-tumor datasets. Many of the benchmarks in our manuscript were done on non-tumor datasets, such PBMC and brain cells. Simply set key=NULL when constructing the prism object. Please refer to the vignette for details.

# Doubts regarding cell state and cell label; log transformed input #15
# https://github.com/Danko-Lab/BayesPrism/issues/15

# Thank you for your interest in our work. To answer your question: cell type is your desired granularity, while cell state is the very fine granularity where there is a great extent of uncertainty that we want to marginalize (for example, to be used to model subclusters of cells lie on a continuum or have some degree of heterogeneity such as the cancer cells). The definition can be a bit problem-specific.
# 
# Please note that cell state is never a required input, if cell state is not provided, BayesPrism will treat cell state one-t-one correspond to cell types.

# cell state label is mandatory ?
# https://github.com/Danko-Lab/BayesPrism/issues/66
# Hi, I would like to use BayesPrism on rnaseq data, is it mandatory to have also a file which specify the cell state label of my single cell data? Or can it work with just the cell type label ?
  
# Yes. To only use cell type, simply specify cell.state.labels = the same input you provided for cell.type.labels.

# Theoretically speaking, the definition of cell states were developed to account for the uncertainty in the reference, particularly for cell type whose distribution is heterogeneous and lies on a continuum,  and hence was not developed for accurate quantification. The reason for this is that suppose cell state B is a linear combination of cell state A and C, even though they are presented in the same amount in the mixture, the inference can be unstable (prone to the data quality of the cells that constitute the cell state), and the sparse prior may also penalize cellstate B to close to 0.

# question about input data in BayesPrim #9
# https://github.com/Danko-Lab/BayesPrism/issues/9

# Hello. Thank you for the amazing tool BayesPrim. I am wondering if I have bulk RNA seq data under A, B, C, and D conditions. But I only have the single cell sequence data under condition A. Can I still feed the BayesPrim with all the data to infer the cell type composition in bulk seq data under condition B, C, and D?
# Although BayesPrism is robust to biological variation in the context of tumor heterogeneity, it is hard to predict the extent to which your experimental condition will affect the performance of BayesPrism. This can be highly problem-specific, which ultimately boils down to the question of how much information from scRNA-seq in condition A is retained in other conditions. That being said, I would recommend the use of highly cell-type specific markers (computed from the scRNA-seq of condition A), by implicitly assuming that these marker genes retain their cell type-specificity in other conditions.

# use sparse matrices instead of full matrices ? #58
# https://github.com/Danko-Lab/BayesPrism/issues/58
# Sorry for the delay. The new version v2.2 now supports input of dgMatrix as the input of scRNA-seq reference.
  
