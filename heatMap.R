#!/usr/bin/Rscript
# Script to compute heat map

## Retrieve expression of DEG
DEG_20 <- rbind(Up[1:10, ], Down[1:10, ])

degExpression <- expr0[DEG_20$Ensembl_ID, ]

degExpression <- tibble::column_to_rownames(degExpression, var = "Ensembl_ID")

# Define row label
rowannot <- cbind(rownames(degExpression), rep(c("Up-regulated","Down-regulated"), each = 10,10))

colnames(rowannot) <- c("Gene", "DEG_direction")


Condition <- rep(c("Data1","Data2"),c(x1, x2))

# Example of Status
# pheno == data frame generated in Microarray_preprocessing.R
ref <- pheno$Status

Phenotype <- ref[match(colnames(degExpression), pheno$sample)]


library(ComplexHeatmap)
library(circlize)

## heat map annotation

# Define a data frame containing column annotations 
dfannot = data.frame(Condition,Phenotype)

# Define annotation with colour

ha = HeatmapAnnotation(df = dfannot, show_legend = c(TRUE, TRUE),
                       col = list(Condition = c("Data1" = "red","Data2" = "blue"),
                                  Phenotype = c("Patient" = "darkcyan", "Control" = "darkgoldenrod3")))

# Draw heat map

Heatmap(degExpression, name = "log2(expr)", show_row_names = F, show_column_names = F, top_annotation = ha, cluster_rows = T, 
        cluster_columns = T, clustering_distance_columns = "euclidean", col = c("navy", "white", "red"))




