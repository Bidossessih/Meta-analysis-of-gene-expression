#!/usr/bin/Rscript
## Script to preprocess Affymetrix GeneChip microarray
###Required input is:###
#1) files containing raw data (.CEL). Each file contains samples from a single study



## Set the directory of raw datasets

inputDir <- './' # working directory

# List of files. Only raw data were listed using pattern option 
allFiles <- list.files(inputDir, full=TRUE, recursive=TRUE, pattern="CEL.gz")
allFiles <- split(allFiles, basename(dirname(allFiles)))

## Install Oligo package if necessary
## source("https://bioconductor.org/biocLite.R")
## biocLite("oligo") #Install Oligo

library(oligo)
Sys.setenv(R_THREADS = 4) # Set the R_THREADSenvironment variable to enable parallel preprocessing

## loadPP: function to preprocess via RMA

loadPP <- function(fns) {
  theTypes <- sapply(fns, oligo:::getCelChipType, TRUE)
  theGrps <- split(names(theTypes), theTypes)
  lapply(theGrps, function(fn) rma(read.celfiles(fn)))
}

# RMA normalization for each dataset
res <- lapply(allFiles, loadPP)

res <- unlist(res)

#Annotate datasets
annot <- data.frame(series=gsub('(.*)\\.(.*)', '\\1', names(res)),
                    platform=gsub('(.*)\\.(.*)', '\\2', names(res)))

annot$filter <- paste('affy', gsub('-', '_', annot$platform), sep='_')

annot$filter <- gsub('hu6800', 'hugenefl', tolower(annot$filter))
annot$array = paste0(gsub("[-_]", "", tolower(annot$platform)), ".db")

#Load metadata
#pheno is a data frame with two columns. 
#column 1 contains sample names
#column 2 contains groups labels
pheno <- read.table("Tab1.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pheno$Status <- factor(pheno$Status, levels = c("Control", "Disease")) # example of studies comparing disease samples to control group
                                                                        
#limma is used to fit the linear model
#Coefficients will be used to compute the correlation matrix of the datasets

library(limma) 

## myf with arguments 'eset' e 'si' (sample information)
## fit linear model using moderate t-test
## return only coefficient per 'Status'
myf <- function(eset, si) {
  sns <- sampleNames(eset)
  intPheno <- si[match(sns, si$sample), ]
  design <- model.matrix(~Status, intPheno)
  fit <- eBayes(lmFit(eset, design))
  #return coefficient
  coef(fit)[, 2]
}

#Apply myf to the normalized results 'res'
limmaRes <- lapply(res, myf, pheno)

array = unique(annot$array)


## install 'array' if necessary
## biocLite(array)
sapply(array, require, character.only=TRUE)

out <- outCoefs <- vector("list", nrow(annot))

#function to compute rowMean of probesets
rowMeanGrp <- function(mat, grp) {
  ones <- matrix(1, nrow = nrow(mat), ncol = 1)
  rowsum(mat, grp)/as.integer(rowsum(ones, grp))
}

# coefficient from limma and Expression matrix
for (i in 1:nrow(annot)) {
  inMat <- exprs(res[[i]]) #Retrieve expression data from eSets
  ps <- rownames(inMat)
  tmp <- toTable(get(gsub("\\.db", "ENSEMBL", annot[i, 'array'])))
  idx <- match(tmp[[1]], ps)
  inMat <- inMat[idx, ]
  out[[i]] <- rowMeanGrp(inMat, tmp[[2]])
  outCoefs[[i]] <- tapply(limmaRes[[i]][idx], tmp[[2]], mean)
  rm(inMat, ps, tmp, idx)
}


#determine common genes if differents microarray platforms were analyzed
#to ensure that variables are comparable
commonPS <- sort(Reduce(intersect, lapply(out, rownames)))

#Retrieve expression of genes common to the platforms
expr0 <- lapply(out, "[", i = commonPS, j = ) 
#sapply(expr0, dim)
expr0 <- do.call(cbind, expr0)

#Retrieve coefficient of common genes
coef0 <- do.call(cbind, lapply(outCoefs, "[", i = commonPS))
colnames(coef0) <- gsub("(.*)\\..*", "\\1", names(res))
#dim(coef0)

#coef0[1:10, 1:6]

#Compute the correlation between datasets
#A quick way to assess the comparability of the datasets

panel.cor <- function(x, y, digits = 4, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  text(0.5, 0.5, txt, cex = 2)
}
panel.sms <- function(...) smoothScatter(..., nrpoints = 0, add = TRUE)

#ylim and xlim can be changed if necessary. It depends on the log expression

pairs(coef0, upper.panel = panel.cor, lower.panel = panel.sms, ylim=c(-2,2), xlim=c(-2,2))

###Meta-Analysis will be performed using the generated expr0
### Save this session to continue the meta-analysis 'meta.R'
### 

