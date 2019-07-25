#!/usr/bin/Rscript
# Script to perform functional gene-set analysis
# Gene set analysis can be performed functional analysis based on
# FAIME (Functional Analysis of Individual Microarray/RNAseq Expression) algrithm
## FAIME improve the interpretability of the gene set analysis
# Heat map is drawn to visualize the results

library("seq2pathway.data")
library("seq2pathway")

# Difine a gene set database

data(MsigDB_C5,package="seq2pathway.data")
class(MsigDB_C5)

# Perform the gene set analysis
dat_gene2path_RNA <- gene2pathway_test(dat=dat1, DataBase=MsigDB_C5,
                                       EmpiricalTest=FALSE, alpha=5, logCheck=FALSE, method="FAIME", na.rm=TRUE)

res <- dat_gene2path_RNA$gene2pathway_result.2 # Alternatively altered pathways can be filtered by the computed p-value

####### Heat map annotation 

Condition <- rep(c("Data1","Data2"),c(x1, x2))

# Example of Status
# pheno == data frame generated in Microarray_preprocessing.R
ref <- pheno$Status

Phenotype <- ref[match(colnames(degExpression), pheno$sample)]


library(ComplexHeatmap)
library(circlize)

## heat map annotation

# Define a data frame containing column annotations 
dfannot <- data.frame(GSEpheno[,2], row.names = GSEpheno$Sample)
rownames(dfannot) = colnames(SCDpathwayDEBP)
colnames(dfannot) <- "Status"
# Define annotation with colour

ha = HeatmapAnnotation(df = dfannot, show_legend = c(TRUE, TRUE),
                       col = list(Condition = c("Data1" = "red","Data2" = "blue"),
                                  Phenotype = c("Patient" = "darkcyan", "Control" = "darkgoldenrod3")))



Heatmap(res, name = "log2(expr)", show_row_names = T, show_column_names = F, cluster_rows = T, 
                   top_annotation = ha, cluster_columns = F, clustering_distance_columns = "euclidean")
  
