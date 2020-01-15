# as for r_frequently_used_rnaseq

library(org.Hs.eg.db)
library("DEGreport")
library("clusterProfiler")
library(pathview)
library(gage)
library(DESeq2)
library("RColorBrewer")
library(EnhancedVolcano)
library(edgeR)
library("gplots")

setwd("~/Documents/...")


# get b2b vs hela:
# load data as before.....
celltype<-factor(c(rep("B2B_Un",3),
                   rep("B2B_T",3),
                   rep("Hela_Un",3),
                   rep("Hela_T",3)))
des<-formula(~celltype)

countsTable <-read.csv("revcounts.csv",check.names=FALSE,row.names=1)
myNames<-colnames(countsTable[c(8:19)])
colDataNames<-data.frame(row.names=myNames, celltype=celltype)
ddsHTSeq<-DESeqDataSetFromMatrix(countsTable[c(8:19)], 
                                 colData=colDataNames, design=des, ignoreRank = FALSE)


#Run DeSeq analysis
dds<-DESeq(ddsHTSeq,betaPrior=FALSE)

sigs <- significants(degs[[2]], fc = 0, fdr = 0.05)
length(sigs)
sigs <- gsub("\\..*", "", sigs)

# subset for conditions we want to compare
untreated <- countsTable[c(8,9,10,14,15,16)]
hela <- countsTable[c(14,15,16,17,18,19)]
b2b <- countsTable[c(8,9,10,11,12,13)]

# get counts per million (from edgeR)
logcounts_untreated <- cpm(untreated,log=TRUE)
logcounts_b2b <- cpm(b2b,log=TRUE)
logcounts_hela <- cpm(hela,log=TRUE)

# make names better for heatmap presentation- remove .x then convert to gene symbol
row.names(logcounts_untreated) <- gsub("\\..*", "", row.names(logcounts_untreated))
row.names(logcounts_untreated) <- mapIds(org.Hs.eg.db, row.names(logcounts_untreated), "SYMBOL", "ENSEMBL")

row.names(logcounts_hela) <- gsub("\\..*", "", row.names(logcounts_hela))
row.names(logcounts_hela) <- mapIds(org.Hs.eg.db, row.names(logcounts_hela), "SYMBOL", "ENSEMBL")

row.names(logcounts_b2b) <- gsub("\\..*", "", row.names(logcounts_b2b))
row.names(logcounts_b2b) <- mapIds(org.Hs.eg.db, row.names(logcounts_b2b), "SYMBOL", "ENSEMBL")

# names of genes of interest
glycocalyx_genes <- c("NANS", "NANP", "CMAS",  "GNE", "NEU4", "NEU1", "NEU2", "NEU3", "SLC35A1", "SLC17A5", "SDC1", "MUC16", "SIGLEC11", "FTL", "FTH1", "MUC1", "CEACAM7", "DEFB126", "MUC7", "HPSE", "MUC4", "ITIH2", "MUC17", "SOD3", "GPC1", "MUC5B", "ZP1", "ZG16", "THBD")
# mapped to ensembl names
glycocalyx_genes_ensembl <- mapIds(org.Hs.eg.db, glycocalyx_genes, "ENSEMBL", "SYMBOL")
sigs_gc <- glycocalyx_genes_ensembl[glycocalyx_genes_ensembl %in% sigs]

# heatmap with rows ordered by similarity to group up/down regulated genes together
heatmap.2(data.matrix(glycocalyx_results),margins=c(10,10), scale="row", labCol=c("B2B_UN_1", "B2B_UN_2", "B2B_UN_3", "HELA_UN_1", "HELA_UN_2", "HELA_UN_3"),col = bluered(100), trace="none", dendrogram="none", reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
  distfun=function(x) as.dist(1-cor(t(x))),
  hclustfun=function(x) hclust(x, method="complete"))

#############  volcano plot


## get res
res <- results(dds, contrast=c("celltype", "B2B_Un", "Hela_Un"),cooksCutoff=FALSE,independentFiltering=FALSE)

row.names(res) <- gsub("\\..*", "", row.names(res))


untreated_results <- res[row.names(res) %in% macropinocytosis_genes_ensembl, ]

rownames(untreated_results) <- mapIds(org.Hs.eg.db, rownames(untreated_results), "SYMBOL", "ENSEMBL")

png("macropinocytosis_untreated_volcano.png")
EnhancedVolcano(untreated_results,
    lab = rownames(untreated_results),
    x = 'log2FoldChange', y='pvalue')

