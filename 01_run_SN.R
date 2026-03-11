library(edgeR)
library(data.table)

source("/mnt/lareaulab/reliscu/code/SampleNetwork/SampleNetwork_1.08.r")

setwd("/mnt/lareaulab/reliscu/projects/NSF_GRFP/analyses/bulk/hahn_2023/cortex")

####################################################################################################################################################
############################################################ Counts data ###########################################################################
####################################################################################################################################################

# Counts data

datExprT <- fread("/mnt/lareaulab/reliscu/projects/NSF_GRFP/data/bulk/hahn_2023/cortex/hahn_2023_cortex_STAR_gene_counts.csv", data.table=FALSE)
sampleinfo1 <- fread("/mnt/lareaulab/reliscu/projects/NSF_GRFP/data/bulk/hahn_2023/cortex/hahn_2023_cortex_sampleinfo.csv", data.table=FALSE)

sampleinfo1$grouplabels1 <- "All"

# Order columns by variable they will be colored by:

sampleinfo1 <- sampleinfo1[order(sampleinfo1$tissueLong),]
datExprT <- datExprT[, c(1, match(sampleinfo1[,1], colnames(datExprT)))]

all.equal(sampleinfo1[,1], colnames(datExprT)[-1])

dim(datExprT)

tmmf <- calcNormFactors(datExprT[,-1], method="TMM")
datExprT[,-1] <- sweep(datExprT[,-1], MARGIN=2, FUN="/", STATS=tmmf)

projectname1 <- "hahn_2023_cortex_STAR_counts_TMMF"
skip1 <- 1 # An integer describing the number of feature information columns
indices1 <- list(seq(2, ncol(datExprT)))
samplelabels1 <- 1 # An integer that points to the column number in sampleinfo1 containing the sample labels that will appear in plots.  Note: these 					sample labels must be identical to the sample column headers in datExprT or an error will be triggered.
grouplabels1 <- grep("grouplabels1", colnames(sampleinfo1))
subgroup1 <- grep("tissueLong", colnames(sampleinfo1)) 

btrait1 <- c(
  grep("age", colnames(sampleinfo1)),
  grep("sex", colnames(sampleinfo1)),
  grep("tissueLong", colnames(sampleinfo1)),
  grep("mouse_id_number", colnames(sampleinfo1))
)

asfactors1 <- c(
  grep("sex", colnames(sampleinfo1)),
  grep("tissueLong", colnames(sampleinfo1)),
  grep("mouse_id_number", colnames(sampleinfo1))
)

SampleNetwork(
  datExprT=datExprT,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=skip1,
  indices1=indices1,
  sampleinfo1=sampleinfo1,
  subgroup1=subgroup1,
  samplelabels1=samplelabels1,
  grouplabels1=grouplabels1,
  fitmodels1=TRUE,
  whichmodel1="univariate",
  whichfit1="pc1",
  btrait1=btrait1,
  trait1=NULL,
  asfactors1=NULL,
  projectname1=projectname1,
  cexlabels1=0.7,
  normalize1=TRUE,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)

####################################################################################################################################################
############################################################### TPM data ###########################################################################
####################################################################################################################################################

datExprT <- fread("/mnt/lareaulab/reliscu/projects/NSF_GRFP/data/bulk/hahn_2023/cortex/hahn_2023_cortex_STAR_TPM_counts.csv", data.table=FALSE)
sampleinfo1 <- fread("/mnt/lareaulab/reliscu/projects/NSF_GRFP/data/bulk/hahn_2023/cortex/hahn_2023_cortex_sampleinfo.csv", data.table=FALSE)

sampleinfo1$grouplabels1 <- "All"

# Order columns by variable they will be colored by:

sampleinfo1 <- sampleinfo1[order(sampleinfo1$tissueLong),]
datExprT <- datExprT[, c(1, match(sampleinfo1[,1], colnames(datExprT)))]

all.equal(sampleinfo1[,1], colnames(datExprT)[-1])

dim(datExprT)

projectname1 <- "hahn_2023_cortex_STAR_TPM"
skip1 <- 1 # An integer describing the number of feature information columns
indices1 <- list(seq(2, ncol(datExprT)))
samplelabels1 <- 1 # An integer that points to the column number in sampleinfo1 containing the sample labels that will appear in plots.  Note: these 					sample labels must be identical to the sample column headers in datExprT or an error will be triggered.
grouplabels1 <- grep("grouplabels1", colnames(sampleinfo1))
subgroup1 <- grep("tissueLong", colnames(sampleinfo1)) 

btrait1 <- c(
  grep("age", colnames(sampleinfo1)),
  grep("sex", colnames(sampleinfo1)),
  grep("tissueLong", colnames(sampleinfo1)),
  grep("mouse_id_number", colnames(sampleinfo1))
)

asfactors1 <- c(
  grep("sex", colnames(sampleinfo1)),
  grep("tissueLong", colnames(sampleinfo1)),
  grep("mouse_id_number", colnames(sampleinfo1))
)

SampleNetwork(
  datExprT=datExprT,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=skip1,
  indices1=indices1,
  sampleinfo1=sampleinfo1,
  subgroup1=subgroup1,
  samplelabels1=samplelabels1,
  grouplabels1=grouplabels1,
  fitmodels1=TRUE,
  whichmodel1="univariate",
  whichfit1="pc1",
  btrait1=btrait1,
  trait1=NULL,
  asfactors1=NULL,
  projectname1=projectname1,
  cexlabels1=0.7,
  normalize1=TRUE,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)