library("ggplot2") # basic plotting
library("RColorBrewer") # colour schemes for plots
library("gplots") # heatmap plotting
library("ggdendro") # dendrogram
library("ggrepel")
library("DESeq2") # differential expression analysis
library("readxl") 
library("genefilter")
library("pcaExplorer") # pca analysis
library(org.Hs.eg.db) # used for converting gene names, e.g. from ensembl to gene symbol
library("DEGreport") # enhanced analysis and reporting of differential expression
library("clusterProfiler") # more gene analysis
library(pathview) # view pathways analysis
library(gage)
library(EnhancedVolcano) # plot volcano 
library(edgeR) # more differential expression analysis


setwd("~/Documents/....")


#### PCA and DESeq lifted from Rachel Scholey's script

# expriment design- example, we have three replicates each of 4 conditions:
# B2B treated and untreated, and HeLa treated and untreated
celltype<-factor(c(rep("B2B_Un",3),
                   rep("B2B_T",3),
                   rep("Hela_Un",3),
                   rep("Hela_T",3)))
des<-formula(~celltype)


# load data from CSV
countsTable <-read.csv("revcounts_only.csv",check.names=FALSE,row.names=1)
# specify count columns for each replicate, rows are gene names
myNames<-colnames(countsTable[c(8:19)])
colDataNames<-data.frame(row.names=myNames, celltype=celltype)
# get data format required for running DESeq
ddsHTSeq<-DESeqDataSetFromMatrix(countsTable[c(8:19)], 
                                 colData=colDataNames, design=des, ignoreRank = FALSE)

#Run DeSeq analysis
dds<-DESeq(ddsHTSeq,betaPrior=FALSE)
resultsNames(dds)
dds$Condition <- celltype
## get normalised counts
normCounts<-as.data.frame(counts(dds,normalized=TRUE))
write.csv(normCounts,"SV_DESeq_Normalised.csv") 

## get "mean" counts per group
Means <- sapply( levels(dds$Condition), 
                     function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$Condition == lvl,drop = FALSE] ) )


dds$Condition <- celltype
write.csv(Means,"SV_DESeq_Means.csv")


## do rlog transformation for PCA analysis
rld <- rlogTransformation(dds, blind=TRUE)

#scree plot
## calculate the variance for each gene
rv <- rowVars(assay(rld))
## select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
## perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(rld)[select,]))
## the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
##plot the "percentVar"
scree_plot=data.frame(percentVar)
scree_plot[,2]<- c(1:12)
colnames(scree_plot)<-c("variance","component_number")
(ggplot(scree_plot, mapping=aes(x=component_number, y=variance))+geom_bar(stat="identity")+
    ggtitle("Scree Plot: Variance Explained by Principle Components"))


#generate PCA plot based on rlogtransformation

###PCA PLOT full dataset PC1 vs 2
PCAdata<-plotPCA(rld, intgroup=c("Condition"),ntop = 500,returnData=TRUE)
##add extra info to PCA data to help annotate PCA plot
PCAdata$CellLine <- c(rep("B2B",6),rep("HeLa",6))
PCAdata$GO_Treatment <- c(rep("Untreated",3), rep("Treated",3), rep("Untreated",3), rep("Treated",3))
PCAdata$sample <- c("1","2","3","4","5","6","7","8","9","10","11","12")
PCAdata # check PCA data table looks ok

## plot PCA plot PC1 vs 2
percentVar <- round(100 * attr(PCAdata, "percentVar"))

PCAplot<-qplot(PC1, PC2, colour=GO_Treatment, shape=CellLine, data=PCAdata) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ggtitle("PCA: PC1 vs PC2") +
  scale_size_manual(values=c(3,5,7,9))+
  geom_point(size=6) +
  geom_text_repel(label=PCAdata$sample) +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  guides(shape = guide_legend(override.aes = list(size=6))) +
  theme(plot.title = element_text(colour="black", size = 16)) + 
  theme(axis.title = element_text(colour="black", size = 14)) + 
  theme(legend.title = element_text(colour="black", size = 16)) + 
  theme(legend.text = element_text(colour="black", size = 16)) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
PCAplot

###############################################################

# to split above for specific subsets
# specify new design
treatment<-factor(c(rep("Un",3),rep("T",3)),levels=c("Un","T"))
des<-formula(~treatment)
# take counts for specific columns of interest - this is b2b
myNames<-colnames(countsTable[c(8:13)])

# now repeat everything above, but with appropriately subsetted counts table
# can do for hela un/t, then b2b un vs hela un etc



