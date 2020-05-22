library(enrichR)

gene_names <- read.csv("gene_ids.txt")

databases<-listEnrichrDbs()
databases<-databases[order(databases$libraryName),3]
database<-c(databases[c(164,154,128,90,68,63,58, 115, 83)])
enrichment<-enrichr(as.character(gene_names$Name), databases = database)

# or $any from database

# https://rdrr.io/github/guokai8/EnrichR/src/R/enrich.R
enrichdot<-function(resultFis,top=50,pvalue.cutoff=0.05,order=FALSE,
                    low="lightpink",high="red",alpha=0.7,
                    padj.cutoff=NULL,fontsize.x=10,fontsize.y=10,
                    usePadj=TRUE,filename=NULL,width=10,height=8){
  library(ggplot2)
  if(!is.null(padj.cutoff)){
    resultFis<-resultFis[resultFis$Padj<padj.cutoff,]
  }else{
    resultFis<-resultFis[resultFis$Pvalue<pvalue.cutoff,]
  }
  if(nrow(resultFis)>=top){
    dd<-resultFis[1:top,]
  }else{
    dd<-resultFis
  }
  if(nrow(dd)>=1){
    if(order==TRUE){
      dd$Term<-factor(dd$Term,levels=dd$Term[order(dd$Overlap)])
    }
    if(usePadj==FALSE){
      p<-ggplot(dd,aes(x=GeneRatio,y=Term))+geom_point(aes(size=Overlap,color=-log10(Pvalue)),alpha=alpha)+
        theme(axis.text.y=element_text(face="bold",size=fontsize.y),axis.text.x=element_text(face="bold",color="black",size=fontsize.x))+
        scale_colour_gradient(low=low,high=high)+theme_minimal()+ylab("Pathway name")+
        xlab("Overlap")+labs(size="Gene number")+guides(color=guide_colourbar(order = 1),size=guide_legend(order = 2))
      print(p)
    }else{
      p<-ggplot(dd,aes(x=GeneRatio,y=Term))+geom_point(aes(size=Overlap,color=-log10(Padj)),alpha=alpha)+
        theme(axis.text.y=element_text(face="bold",size=fontsize.y),axis.text.x=element_text(face="bold",color="black",size=fontsize.x))+
        scale_colour_gradient(low=low,high=high)+theme_minimal()+ylab("Pathway name")+
        xlab("Overlap factor")+labs(size="Gene number")+guides(color=guide_colourbar(order = 1),size=guide_legend(order = 2))
      print(p)
    }
      if(!is.null(filename)){
        ggsave(p,file=paste(filename,"dotplot.pdf",sep="_"),width=width,height=height)
      }
    }else{
      cat("No Pathway enrichment results were found!\n")
    }

}


# or $any from database
df <- enrichment$Jensen_DISEASES
# calculate x axis value
df$GeneRatio <- sapply(df$Overlap, function(x) eval(parse(text=x)))
names(df)[names(df) == "Adjusted.P.value"] <- "Padj"
names(df)[names(df) == "P.value"] <- "Pvalue"

enrichdot(df, padj.cutoff=0.05, usePadj=TRUE, order=TRUE, top=15)

