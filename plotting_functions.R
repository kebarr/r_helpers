library("ggplot2") # basic plotting
library("RColorBrewer") # colour schemes for plots
library("gplots") # heatmap plotting
library("ggdendro") # dendrogram
library("ggrepel")
library("DESeq2") # differential expression analysis
library("readxl") 
library("genefilter")
library("pcaExplorer") # pca analysis


#### LIFTED FROM RACHELS SCRIPT, NOT MY WORK!!!!


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#############################################
# helper function for creating dendograms
#https://plot.ly/ggplot2/ggdendro-dendrograms/
############################################
ggdend <- function(df) {
  ggplot() +
    geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
    labs(x = "", y = "") + theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())
}
#dendogram
#x<- t(assay(vsd)[select,])
#dd.col <- as.dendrogram(hclust(dist(x)))
#dd.row <- as.dendrogram(hclust(dist(t(x))))
#dx <- dendro_data(dd.row)
#dy <- dendro_data(dd.col)
# x/y dendograms
#px <- ggdend(dx$segments)
#myDendro_y_PC <- ggdend(dy$segments) + coord_flip()

####################################
#function my code edit of plotPCA
####################################
plotPCALeo<-function (x, intgroup = "Treatment", ntop = 500, returnData = FALSE, PCx=1, PCy=2)
{
  #rv <- rowVars(assay(x))
  rv = apply((assay(x)), 1, var)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(x)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(x)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(x)[, intgroup, drop = FALSE])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  d <- data.frame(PCX = pca$x[, PCx], PCY = pca$x[, PCy], group = group, 
                  intgroup.df, names = colnames(x))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[PCx:PCy]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PCX", y = "PCY", color = "group")) + 
    #ggplot(data = d, aes_string(x = "PCX", y = "PCY", color=Tgfb1, shape=Treatment)) + 
    geom_point(size = 3) + xlab(paste0("PC",PCx,": ", round(percentVar[1] * 
                                                              100), "% variance")) + ylab(paste0("PC",PCy,": ", round(percentVar[2] * 
                                                                                                                        100), "% variance"))
}
