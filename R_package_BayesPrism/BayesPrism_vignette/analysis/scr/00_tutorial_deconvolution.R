# AIM ---------------------------------------------------------------------
# run the tutorial analysis for the BayesPrism package
# usethis::create_github_token()
# gitcreds::gitcreds_set()

# libraries ---------------------------------------------------------------
library(tidyverse)
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

# bk.dat: The sample-by-gene raw count matrix of bulk RNA-seq expression. rownames are bulk sample IDs, while colnames are gene names/IDs.
dim(bk.dat)
head(rownames(bk.dat))
head(colnames(bk.dat))

# sc.dat: The cell-by-gene raw count matrix of single cell RNA-seq expression. rownames are cell IDs, while colnames are gene names/IDs. sc.dat should be a dense matrix.
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

table(cbind.data.frame(cell.state.labels, cell.type.labels))

# Please make sure that all cell states contain a reasonable number of cells, e.g. >20 or >50, so that their profile can be represented accurately.

# What to supply for cell.state.labels and cell.type.labels? The definition of cell type and cell state can be somewhat arbitrary (similar to the issue of assigning cell types for scRNA-seq) and depends on the question of interest. Their definitions depend on the granularity we aim at and the confidence of the cell.type.labels in scRNA-seq data. Usually, a good rule of thumb is as follows.
# 1) Define cell types as the cluster of cells having a sufficient number of significantly differentially expressed genes than other cell types, e.g., greater than 50 or even 100. For clusters that are too similar in transcription, we recommend treating them as cell states, which will be summed up before the final Gibbs sampling. Therefore, cell states are often suitable for cells that form a continuum on the phenotypic manifold rather than distinct clusters.
# 2) Define multiple cell states for cell types of significant heterogeneity, such as malignant cells, and of interest to deconvolve their transcription.


# QC of cell type and state labels ----------------------------------------
# We recommend first plotting the pairwise correlation matrix between cell states and between cell types. This will give us a sense of their quality. In cases where cell types/states are not represented by sufficient amount of information (low cell count and/or low library size), the low-quality cell types/states tend to cluster together. Users may re-cluster the data at higher granularity, or merge those cell types/states with the most similar cell types/states, or remove them (if re-clustering and merging is not appropriate).

plot.cor.phi (input = sc.dat,
              input.labels = cell.state.labels,
              title="cell state correlation",
              #specify pdf.prefix if need to output to pdf
              #pdf.prefix="gbm.cor.cs", 
              cexRow=0.2, cexCol=0.2,
              margins=c(2,2))


plot.cor.phi (input = sc.dat, 
              input.labels = cell.type.labels, 
              title="cell type correlation",
              #specify pdf.prefix if need to output to pdf
              #pdf.prefix="gbm.cor.ct",
              cexRow=0.5, cexCol=0.5,
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

sc.dat.filtered <- cleanup.genes (input = sc.dat,
                                  input.type = "count.matrix",
                                  species = "hs", 
                                  gene.group = c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells = 5)

dim(sc.dat.filtered)

# Next, we check the concordance of gene expression for different types of genes. As bulk and single cell data are usually collected by different experimental protocols, they may have different sensitivity to different types of genes.
# note this function only works for human data. For other species, you are advised to make plots by yourself.
plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = bk.dat
                 #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
                 )

# We observe that protein coding genes are the most concordant group between two assays. To reduce batch effects and speed up computation, we perform deconvolution on protein coding genes. We have also tried to run BayesPrism on all genes. The results were similar.

# Subset protein coding genes.
sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
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

sc.dat.filtered.pc.sig <- select.marker (sc.dat = sc.dat.filtered.pc,
                                         stat = diff.exp.stat,
                                         pval.max = 0.01,
                                         lfc.min = 0.1)

dim(sc.dat.filtered.pc.sig)

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

Report the summary statistics.

bp.res
#> Input prism info: 
#> Cell states in each cell type: 
#> $tumor
#>  [1] "PJ016-tumor-0" "PJ016-tumor-3" "PJ016-tumor-2" "PJ016-tumor-6" "PJ016-tumor-4" "PJ016-tumor-1" "PJ016-tumor-5" "PJ017-tumor-3" "PJ017-tumor-2" "PJ017-tumor-1" "PJ017-tumor-0" "PJ017-tumor-4" "PJ017-tumor-5" "PJ017-tumor-6" "PJ018-tumor-4" "PJ018-tumor-2" "PJ018-tumor-3" "PJ018-tumor-1" "PJ018-tumor-0" "PJ018-tumor-5" "PJ025-tumor-9" "PJ025-tumor-2" "PJ025-tumor-4" "PJ025-tumor-3" "PJ025-tumor-0" "PJ025-tumor-1" "PJ025-tumor-8" "PJ025-tumor-7" "PJ025-tumor-5" "PJ025-tumor-6" "PJ030-tumor-4" "PJ030-tumor-0" "PJ030-tumor-1" "PJ030-tumor-5" "PJ030-tumor-3" "PJ030-tumor-2" "PJ032-tumor-0" "PJ032-tumor-1" "PJ032-tumor-5" "PJ032-tumor-2" "PJ032-tumor-4" "PJ032-tumor-3" "PJ035-tumor-2" "PJ035-tumor-3" "PJ035-tumor-6" "PJ035-tumor-1" "PJ035-tumor-0" "PJ035-tumor-4" "PJ035-tumor-5" "PJ035-tumor-8" "PJ035-tumor-7" "PJ048-tumor-0" "PJ048-tumor-3" "PJ048-tumor-4" "PJ048-tumor-2" "PJ048-tumor-6" "PJ048-tumor-1" "PJ048-tumor-8" "PJ048-tumor-7" "PJ048-tumor-5"
#> 
#> $myeloid
#> [1] "myeloid_2" "myeloid_5" "myeloid_3" "myeloid_7" "myeloid_0" "myeloid_6" "myeloid_4" "myeloid_1" "myeloid_8"
#> 
#> $pericyte
#> [1] "pericyte"
#> 
#> $endothelial
#> [1] "endothelial"
#> 
#> $tcell
#> [1] "tcell"
#> 
#> $oligo
#> [1] "oligo"
#> 
#> 
#> Identifier of the malignant cell type:  tumor 
#> Number of cell states:  73 
#> Number of cell types:  6 
#> Number of mixtures:  169 
#> Number of genes:  16145 
#> 
#> Initial cell type fractions: 
#>         tumor myeloid pericyte endothelial tcell oligo
#> Min.    0.199   0.004    0.000       0.015  0.00 0.000
#> 1st Qu. 0.702   0.051    0.014       0.036  0.00 0.010
#> Median  0.775   0.102    0.025       0.046  0.00 0.025
#> Mean    0.758   0.109    0.043       0.048  0.00 0.041
#> 3rd Qu. 0.848   0.145    0.048       0.060  0.00 0.050
#> Max.    0.952   0.578    0.503       0.133  0.01 0.329
#> Updated cell type fractions: 
#>         tumor myeloid pericyte endothelial tcell oligo
#> Min.    0.272   0.000    0.000       0.008 0.000 0.000
#> 1st Qu. 0.751   0.046    0.006       0.027 0.000 0.002
#> Median  0.815   0.090    0.014       0.039 0.000 0.011
#> Mean    0.801   0.098    0.030       0.041 0.001 0.029
#> 3rd Qu. 0.878   0.130    0.033       0.053 0.001 0.036
#> Max.    0.968   0.565    0.385       0.134 0.010 0.307
Alternatively, if user would like to run step 4 separately from 1-3, one may do the following:
  
  # not run
  bp.res.initial <- run.prism(prism = myPrism, n.cores=50, update.gibbs=FALSE)

bp.res.update <- update.theta (bp = bp.res.initial)
Extract results
Understand the results BayesPrism is expected to generate the following results: 
  The output from bp.res is an S4 object of the class "BayesPrism". Let’s take a look at what’s inside:
  
  slotNames(bp.res)
#> [1] "prism"                       "posterior.initial.cellState" "posterior.initial.cellType"  "reference.update"            "posterior.theta_f"           "control_param"
"prism" is the input prism.
"posterior.initial.cellState" is the result of step 2.
"posterior.initial.cellType" is the result of step 3.
"reference.update" is the updated reference ψ
.
"posterior.theta_f" is the result of step 4.
"control_param" contains the parameters to run BayesPrism.
We provide utility functions to extract deconvolved cell type fraction and expression from the output of run.prism.

# extract posterior mean of cell type fraction theta
theta <- get.fraction (bp=bp.res,
                       which.theta="final",
                       state.or.type="type")

head(theta)
#>                                  tumor    myeloid     pericyte endothelial        tcell        oligo
#> TCGA-06-2563-01A-01R-1849-01 0.8392297 0.04329259 2.999022e-02  0.07528272 6.474488e-04 0.0115573149
#> TCGA-06-0749-01A-01R-1849-01 0.7090654 0.17001073 8.995526e-07  0.01275526 1.179331e-06 0.1081665709
#> TCGA-06-5418-01A-01R-1849-01 0.8625322 0.09839143 9.729268e-03  0.02416954 7.039913e-07 0.0051768589
#> TCGA-06-0211-01B-01R-1849-01 0.8893449 0.04482991 1.131622e-02  0.05435490 2.508238e-06 0.0001515524
#> TCGA-19-2625-01A-01R-1850-01 0.9406438 0.03546026 1.932740e-03  0.01309753 4.535897e-06 0.0088610997
#> TCGA-19-4065-02A-11R-2005-01 0.6763166 0.08374439 1.849921e-01  0.01918132 3.541126e-04 0.0354114398
# extract coefficient of variation (CV) of cell type fraction
theta.cv <- bp.res@posterior.theta_f@theta.cv

head(theta.cv)
#>                                     tumor      myeloid     pericyte endothelial      tcell       oligo
#> TCGA-06-2563-01A-01R-1849-01 0.0001722829 0.0016025331 0.0026568433 0.001853843 0.05452368 0.005606057
#> TCGA-06-0749-01A-01R-1849-01 0.0002333853 0.0006859332 0.8109050326 0.005683516 0.74729658 0.001276784
#> TCGA-06-5418-01A-01R-1849-01 0.0001601128 0.0009389782 0.0070872942 0.004131925 0.83670176 0.011362855
#> TCGA-06-0211-01B-01R-1849-01 0.0001175529 0.0014412303 0.0055064706 0.001673105 0.60905434 0.280609474
#> TCGA-19-2625-01A-01R-1850-01 0.0001327408 0.0023661566 0.0335184780 0.006286823 0.73453327 0.006138184
#> TCGA-19-4065-02A-11R-2005-01 0.0002862600 0.0012227981 0.0009100677 0.005660069 0.06371147 0.002243368
theta.cv quantifies how concentrated the posterior distribution is. Higher theta value is associated with lower CV. When performing rank-ordered statistics, such as Spearman’s rank correlation, users may mask theta (e.g. clip at zero for small theta) where the CV is above some number, say > 0.1. Lower sequencing depth, which is the case for spatial transcriptomics, is generally associated with higher CV. A reasonable cutoff of CV for deconvolving Visium data may be ~0.5.
# extract posterior mean of cell type-specific gene expression count matrix Z  
Z.tumor <- get.exp (bp=bp.res,
                    state.or.type="type",
                    cell.name="tumor")

head(t(Z.tumor[1:5,]))
#>                    TCGA-06-2563-01A-01R-1849-01 TCGA-06-0749-01A-01R-1849-01 TCGA-06-5418-01A-01R-1849-01 TCGA-06-0211-01B-01R-1849-01 TCGA-19-2625-01A-01R-1850-01
#> ENSG00000130876.10                       55.980                      444.980                        9.000                       56.996                       15.996
#> ENSG00000165244.6                       206.344                       38.096                      291.872                      530.732                      507.788
#> ENSG00000173597.7                        17.100                        6.588                       21.804                       28.744                       10.800
#> ENSG00000158022.6                         0.000                        0.476                        2.616                        0.960                        7.824
#> ENSG00000167220.10                     2273.824                     1060.976                     2775.180                     2368.016                     2084.620
#> ENSG00000126106.12                      678.468                      685.948                      481.000                      698.888                     1296.080
#save the result
save(bp.res, file="bp.res.rdata")
Downstream analysis
Potential downstream analysis can be performed using theta and Z are:
  
  Clustering bulk samples by theta or Z (Z can be normalized by vst(round(t(Z.tumor))), using the vst function from the DESeq2 package.)
Computing z-scores of signature genes for Z of the cell type of interest.
Survival analysis and correlation with other clinical covariates using theta and Z.
Correlating Z (after normalization using vst or from bp.res@reference.update@psi_mal) with theta to understand how gene expression of each gene (in malignant cells) correlates with the cell type fraction of non-malignant cells in tumor micro-environment, followed by gene set enrichment analysis (as done in BayesPrism paper).
Embedding learning of malignant gene expression (see “Tutorial: embedding learning of malignant cell expression using BayesPrism”)