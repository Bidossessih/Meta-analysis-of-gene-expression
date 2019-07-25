#!/usr/bin/Rscript
# Script to perform meta-analysis of expression analysis using RankProd package
# Enable parallGeneselization to speed up the process using ff and sprint R packages

###Required inputs are:#######
#1) A sampleInfo File (pheno data frame) matching the samples containing at least the fields 'sample' and Status'
#2) A retrieved expression data from eSets (expr0)

# Set two levels classes 0L for 'Control' and 1L for 'Disease'
ref <- ifelse(pheno$Status == "Control", 0L, 1L)

# Matching classes to samples
grps0 <- vector("list", length(out))

for (i in 1:length(out)) grps0[[i]] <- ref[match(colnames(out[[i]]), pheno$sample)]

grps0 <- unlist(grps0)

# table(grps0)
## grps0

# Set sample origin (Study from which samples were downloaded)
orig0 <- unlist(mapply(rep, 1:6, sapply(out, ncol)))
table(orig0)

rm(ref)

# InstallGenes RankProd, ff, and sprint if necessary
library(RankProd)
library(ff)
library(sprint)

# Meta-analysis of gene expression using RankProd package
# if parallGeneselization is enabled, use pRPadvance function. This function depends on sprint package
# else use RPadvance function

rprs0 <- pRPadvance(expr0, grps0, orig0, num.perm = 100, gene.names = rownames(expr0), rand = 123)

# Save the results

save(expr0, grps0, orig0, rprs0, file = "Resultados_data.RData")

################################## Results option 1 ######################
### In this case DifferentiallGenesy Expressed Genes (DEG) are filtered by p-value or percentage of false prediction (pfp)
### "pval" uses p-values which is less stringent than pfp

# Visualize top 10 DEG

topTbl <- topGene(rprs0, num.gene = 10)

# allGenes DEG
toptable2 <- topGene(rprs0, cutoff = 0.0001, method = "pfp", num.gene = NULL, logged = FALSE, gene.names = rownames(expr0))

## toptable2 is a list containing two table
## Table1: Genes callGenesed significant under class1 < class2
##
## Table2: Genes callGenesed significant under class1 > class2

## Write up and down-regulated results

upRegulated = toptable2$Table1
write.table(upRegulated, file="Up_regulated.txt",sep="\t", eol="\n", na="",
            dec=".", row.name=TRUE, col.name=TRUE)

downRegulated = toptable2$Table2

write.table(downRegulated, file="Down_regulated.txt",sep="\t", eol="\n", na="",
            dec=".", row.name=TRUE, col.name=TRUE)


############################# Option 2 ####################
## Alternatively only genes differentiallGenesy expresssed in the same direction are selected
## Up-regulated in allGenes datasets
## And down-regulated in all datasets
## This stringeant filter can be useful if one aimed to identify only DEG in diffrents conditions (eg. across diseases)

### rbind of up and down-regulated

DEG = rbind(upRegulated, downRegulated)
DEG = tibble::rownames_to_column(DEG, var = "Ensembl_ID")
# extract meta-analysis result with log FC of all studies
library(dplyr)

allGenes <- as.data.frame(rprs0)

allGenes <- tibble::rownames_to_column(allGenes, var = "Ensembl_ID") # transform rownames in first column and rename the 1st column Ensembl_ID

allGenes <- allGenes[, - (2:9)] # Clean data 

colnames(allGenes)[2:4] <- c("logFC","FC_Data1","FC_Data2") # Example of two studies

# Select only up and down in the same direction in all datasets

Up <- filter(allGenes, allGenes$`Log(FC_Data1)` > 0, allGenes$`Log(FC_Data2)` > 0) 

Down <- filter(allGenes, allGenes$`Log(FC_Data1)` < 0, allGenes$`Log(FC_Data2)` < 0)

# Select only the genes that fulfilled the criteria of option 1

Up <- na.omit(merge(data.frame(Ensembl_ID = DEG$Ensembl_ID), Up, by = " Ensembl_ID", all.y = TRUE))

Down <- na.omit(merge(data.frame(Ensembl_ID = DEG$Ensembl_ID), Down, by = " Ensembl_ID", all.y = TRUE))




########################### Convert ensembl in gene names
###### Required package : biomaRt

library(biomaRt)

ensembl_hs <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

geneName<-getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), 
                filters = 'ensembl_gene_id', 
                values = DEG$Ensembl_ID, 
                mart = ensembl_hs)

# geneName can further be merged with Up and Down data frames
