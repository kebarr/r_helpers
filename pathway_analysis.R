library(org.Hs.eg.db)
library("DEGreport")
library("clusterProfiler")
library(pathview)
library(gage)
library(DESeq2)




# set up standard design- in this case hela treated vs untreated
treatment<-factor(c(rep("Un",3),rep("T",3)),levels=c("Un","T"))
des<-formula(~treatment)
countsTable <-read.csv("revcounts_only.csv",check.names=FALSE,row.names=1)

# relevant columns from csv
myNames_hela<-colnames(countsTable[c(14:19)]) 
colDataNames_hela<-data.frame(row.names=myNames_hela, treatment=treatment)
ddsHTSeq_hela<-DESeqDataSetFromMatrix(countsTable[c(14:19)], 
                                 colData=colDataNames_hela, design=des, ignoreRank = FALSE)

#Run DeSeq analysis
dds_hela<-DESeq(ddsHTSeq_hela,betaPrior=FALSE)
resultsNames(dds_hela)

# some info about this function here: https://bioconductor.riken.jp/packages/devel/bioc/vignettes/DEGreport/inst/doc/DEGreport.html

degs_hela <- degComps(dds_hela, combs = "treatment",
                 contrast = list("treatment_T_vs_Un",
                                 c("treatment", "T", "Un")))


# get significantly differentially expressed genes
sigs_hela <- significants(degs_hela[[2]], fc = 0, fdr = 0.05) # one and 2 are same for data like this
length(sigs_hela) 
sigs_hela_gn <- gsub("\\..*", "", sigs_hela)





