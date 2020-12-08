#source("https://bioconductor.org/biocLite.R")
#.libPaths( "C:/Program Files/R/R-3.3.3/library")
library(fastmatch)
library(reshape2)
#source("R/general_utils.R")
#options(java.parameters = "-Xmx8000m")
#library(httRws)
library(openxlsx)
library(tidyr)
library(grDevices)
library(RColorBrewer)
library(gplots)
library(ape)
#--------------------------------------------------------------------------------------
#' Prepare the raw data for the replicated chemicals
#'
#' @param dir The directory to work in
#'
#--------------------------------------------------------------------------------------
raw.control.analysis.2 <- function(to.file=F,
                                   do.read=F,
                                   do.prep=F,
                                   do.step.1=F,
                                   do.step.2=F,
                                   do.step.3=F,
                                   do.step.4=F,
                                   dirin="../input/httr_mcf7_screen/raw_ctrl/",
                                   welltype="dmso",
                                   well_type_prefix="",
                                   dirout="../output/raw_count_analysis/",
                                   cutoff=5) {
  printCurrentFunction()
  if(do.read) {
    file <- paste0(dirin,"httr-ph-i-raw-pl-bl_",welltype,"-v1.tsv")
    print(file)
    mat <- read.delim(file,stringsAsFactors=F)
    x <- unite(mat[,c("plate_id","block_id","pg_id","well_id")],"rowkey")
    mat$rowkey <- x[,1]
    CTRLMAT <<- mat
    pmat <- unique(CTRLMAT[,c("plate_id","block_id","pg_id","well_id")])
    x <- unite(pmat,"rowkey")
    rownames(pmat) <- x[,1]
    PLATEMAT <<- pmat
  }

  if(do.prep) {
    mat <- CTRLMAT[,c("rowkey","probe_id","probe_count")]
    countmat <- dcast(mat,rowkey~probe_id)
    rownames(countmat) <- countmat[,"rowkey"]
    countmat <- countmat[,2:ncol(countmat)]
    imat <- countmat
    imat[!is.na(imat)] <- 1
    imat[is.na(imat)] <- 0
    cs <- colSums(imat)
    cs <- cs/nrow(countmat)
    countmat2 <- countmat[,cs>0.95]
    countmat2[is.na(countmat2)] <- 0.5
    countmat2 <- log2(countmat2)
    COUNTMAT <<- countmat2
  }
  if(do.step.1) {
    countmat <- as.matrix(COUNTMAT)
    if(to.file) {
      fname <- paste0(dirout,"raw_count_comparison_",welltype,"_step_1.pdf")
      pdf(file=fname,width=7,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    }
    par(mfrow=c(3,2),mar=c(4,4,2,2))
    pg.list <- sort(unique(PLATEMAT$pg_id))
    for(i in 1:length(pg.list)) {
      pg_id <- pg.list[i]
      block <- unique(PLATEMAT[is.element(PLATEMAT[,"pg_id"],pg_id),"block_id"])
      cat("Step 1:",pg_id,":",block,"\n")
      row.list <- rownames(PLATEMAT[is.element(PLATEMAT[,"pg_id"],pg_id),])
      temp <- countmat[row.list,]
      meds <- apply(temp,FUN=median,MARGIN=2)
      mads <- apply(temp,FUN=mad,MARGIN=2)
      plot(mads~meds,xlab="med (log2 count)",ylab="mad (log2 counts)",main=paste("plate group",pg_id," block",block),pch=21,cex=0.2,xlim=c(0,15),ylim=c(0,5))
      if(!to.file) browser()
    }
    if(to.file) dev.off()

  }
  if(do.step.2) {
    countmat <- as.matrix(COUNTMAT)
    if(to.file) {
      fname <- paste0(dirout,"raw_count_comparison_",welltype,"_step_2.pdf")
      pdf(file=fname,width=7,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    }
    par(mfrow=c(4,3),mar=c(2,2,1,1))
    pg.list <- sort(unique(PLATEMAT$pg_id))
    for(i in 1:length(pg.list)) {
      pg_id <- pg.list[i]
      block <- unique(PLATEMAT[is.element(PLATEMAT[,"pg_id"],pg_id),"block_id"])
      cat("Step 2:",pg_id,":",block,"\n")
      row.list <- rownames(PLATEMAT[is.element(PLATEMAT[,"pg_id"],pg_id),])
      nrow <- length(row.list)
      for(j in 1:nrow) {
        names <- NULL
        values <- NULL
        xj <- countmat[row.list[j],]
        for(k in 1:nrow) {
          xk <- countmat[row.list[k],]
          vtemp <- xk-xj
          ntemp <- vtemp
          ntemp[] <- row.list[k]
          names <- c(names,ntemp)
          values <- c(values,vtemp)
        }
        boxplot(values~names,main=row.list[j],ylim=c(-10,10),names=NULL)
        lines(c(-100,100),c(0,0))
        if(!to.file) browser()
      }
       #plot(mads~meds,xlab="med (log2 count)",ylab="mad (log2 counts)",main=paste("plate group",pg_id," block",block),pch=21,cex=0.2,xlim=c(0,15),ylim=c(0,5))
    }
    if(to.file) dev.off()

  }
  if(do.step.3) {
    countmat <- as.matrix(COUNTMAT)
    if(to.file) {
      fname <- paste0(dirout,"raw_count_comparison_",welltype,"_step_3.pdf")
      pdf(file=fname,width=7,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    }

    pmat <- PLATEMAT[rownames(countmat),]
    col.list <- paste0("B",pmat[,"block_id"])
    col.list[is.element(col.list,"B1")] <- "red"
    col.list[is.element(col.list,"B2")] <- "black"
    col.list[is.element(col.list,"B3")] <- "red"
    col.list[is.element(col.list,"B5")] <- "cyan"
    x <- countmat#[,1:1000]
    x[x<0] <- 0
    heatmap.2(x,
              Rowv=T,
              Colv=T,
              hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
              dendrogram="both",
              symm=F,
              scale="none",
              col=brewer.pal(9,"Reds"),
              margins=c(5,8),
              cexCol=0.1,
              cexRow=0.2,
              key=T,
              main="DMSO counts",
              xlab="",
              ylab="",
              trace="none",
              key.title="Key",
              key.xlab="Tlog2(counts)",
              cex.main=1,
              RowSideColors=col.list)
    if(!to.file) browser()
    if(to.file) dev.off()
  }
}
