library(org.Hs.eg.db)
library("DEGreport")
library("clusterProfiler")
library(pathview)
library(gage)
library(DESeq2)
library(pathfindR)

######## GAGE pathway analysis (pathfindR below)
# original gage script is pathways_analysis_checking
# see here:
# https://www.r-bloggers.com/tutorial-rna-seq-differential-expression-pathway-analysis-with-sailfish-deseq2-gage-and-pathview/

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

row.names(degs_hela[[2]]$raw) <- gsub("\\..*", "", row.names(degs_hela[[2]]$raw))
row.names(degs_hela[[2]]$shrunken) <- gsub("\\..*", "", row.names(degs_hela[[2]]$shrunken))

countsTable <-read.csv("revcounts_only.csv",check.names=FALSE,row.names=1)
myNames<-colnames(countsTable[c(8:19)])

# get raw counts
counts_only <- countsTable[myNames]
counts_data <- as.matrix(counts_only)

nonzero_counts = rowSums(counts_data) != 0
counts_data = counts_data[nonzero_counts,]
libsizes=colSums(counts_data)
size.factor=libsizes/exp(mean(log(libsizes)))
counts_data_norm=t(t(counts_data)/size.factor)
cnts.norm=log2(counts_data_norm+8)
row.names(cnts.norm) <- gsub("\\..*", "", row.names(cnts.norm))
row.names(cnts.norm) <- mapIds(org.Hs.eg.db, as.character(row.names(cnts.norm)), "ENTREZID", "ENSEMBL")

data(kegg.gs)
# specify which columns are which
ref_hela=7:9
treat_hela=10:12
counts_kegg_pways_hela <- gage(cnts.norm, gsets = kegg.gs, ref = ref_hela,samp = treat_hela, compare ="unpaired")

# info on using gage now
# https://www.r-bloggers.com/tutorial-rna-seq-differential-expression-pathway-analysis-with-sailfish-deseq2-gage-and-pathview/

# differences between the counts for treatment and control
cnts.d_hela= cnts.norm[, treat_hela]-rowMeans(cnts.norm[, ref_hela])
sel_hela <- counts_kegg_pways_hela$greater[, "q.val"] < 0.1 & !is.na(counts_kegg_pways_hela$greater[,"q.val"])
# most upregulated pathways
path.ids_hela <- rownames(counts_kegg_pways_hela$greater)[sel_hela]
path.ids2_hela <- substr(path.ids_hela, 1, 8)
# visualise 
pv.out.list_hela <- sapply(path.ids2_hela[1:3], function(pid) pathview(
                       gene.data = cnts.d_hela, pathway.id = pid,
                       species = "hsa"))
#down-regulated pathways  (top 3) visualized by pathview
sel.l_hela <- counts_kegg_pways_hela$less[, "q.val"] < 0.1 &
            !is.na(counts_kegg_pways_hela$less[,"q.val"])
path.ids.l_hela <- rownames(counts_kegg_pways_hela$less)[sel.l_hela]
path.ids.l2_hela <- substr(path.ids.l_hela, 1, 8)
#visualise
pv.out.list.l_hela <- sapply(path.ids.l2_hela[1:3], function(pid) pathview(
                    gene.data = cnts.d_hela, pathway.id = pid,
                    species = "hsa"))


######### pathfindR

# some useful info here: https://www.biostars.org/p/322415/


# get input, and add required column names
hela <- degs_hela[[2]]$raw
hela$Gene.symbol <- mapIds(org.Hs.eg.db, as.character(row.names(hela)), "SYMBOL", "ENSEMBL")
hela$logFC <- hela$log2FoldChange
hela$adj.P.Val <- hela$padj

# cols needed by pathfindR
cols <- c("Gene.symbol", "logFC", "adj.P.Val")
hela_final <- na.exclude(data.frame(hela[cols]))

# this takes quite a while.....
hela_sig_output <- run_pathfindR(hela_sig_final)



