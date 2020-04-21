library(fastmatch)
library(reshape2)
#library(httRws)
library(openxlsx)
library(tidyr)
library(grDevices)
library(RColorBrewer)
library(gplots)
library(ape)
#source("R/general_utils.R")
#options(java.parameters = "-Xmx8000m")
#--------------------------------------------------------------------------------------
#' Prepare the raw data for the replicated chemicals
#'
#' @param dir The directory to work in
#'
#--------------------------------------------------------------------------------------
raw.data.prep <- function(welltype="refchems",
                          dirin="../../../Raw/ctrl/",
                          dirout="analysis/raw_count_analysis/data/",
                          cutoff=0.95) {
  printCurrentFunction()
  if(welltype=="refchems") {
    file <- paste0(dirin,"httr-ph-i-raw-pl-",welltype,"-v1.tsv")
    print(file)
    mat.in <- read.delim(file,stringsAsFactors=F)
    chem.list <- c("TSA","GEN","SIRO")
    for(chem in chem.list) {
      print(chem)
      mat <- mat.in[is.element(mat.in[,"name"],chem),]
      x <- unite(mat[,c("plate_id","block_id","pg_id","well_id")],"rowkey")
      mat$rowkey <- x[,1]
      CTRLMAT <- mat
      pmat <- unique(CTRLMAT[,c("plate_id","block_id","pg_id","well_id")])
      x <- unite(pmat,"rowkey")
      rownames(pmat) <- x[,1]
      PLATEMAT <- pmat

      mat <- CTRLMAT[,c("rowkey","probe_id","probe_count")]
       print(dim(mat))
      cat("    start casting\n")
      countmat <- dcast(mat,rowkey~probe_id)
      cat("   casting is done\n")
      rownames(countmat) <- countmat[,"rowkey"]
      countmat <- countmat[,2:ncol(countmat)]
      countmat <- as.matrix(countmat)
      imat <- countmat
      imat[!is.na(imat)] <- 1
      imat[is.na(imat)] <- 0
      cs <- colSums(imat)
      cs <- cs/nrow(countmat)
      countmat <- countmat[,cs>cutoff]
      countmat[is.na(countmat)] <- 0.5
      for(i in 1:nrow(countmat)) countmat[i,] <- 1000000*countmat[i,]/sum(countmat[i,])
      countmat <- log2(countmat)
      COUNTMAT <- countmat
      file <- paste0(dirout,"data_",chem,"_",cutoff,".RData")
      save(CTRLMAT,PLATEMAT,COUNTMAT,file=file)
    }
  }
  else {
    file <- paste0(dirin,"httr-ph-i-raw-pl-",welltype,"-v1.tsv")
    print(file)
    mat <- read.delim(file,stringsAsFactors=F)
    x <- unite(mat[,c("plate_id","block_id","pg_id","well_id")],"rowkey")
    mat$rowkey <- x[,1]
    CTRLMAT <- mat
    pmat <- unique(CTRLMAT[,c("plate_id","block_id","pg_id","well_id")])
    x <- unite(pmat,"rowkey")
    rownames(pmat) <- x[,1]
    PLATEMAT <- pmat

    mat <- CTRLMAT[,c("rowkey","probe_id","probe_count")]
    print(dim(mat))
    cat("    start casting\n")
    countmat <- dcast(mat,rowkey~probe_id)
    cat("   casting is done\n")
    rownames(countmat) <- countmat[,"rowkey"]
    countmat <- countmat[,2:ncol(countmat)]
    countmat <- as.matrix(countmat)
    imat <- countmat
    imat[!is.na(imat)] <- 1
    imat[is.na(imat)] <- 0
    cs <- colSums(imat)
    cs <- cs/nrow(countmat)
    countmat <- countmat[,cs>cutoff]
    countmat[is.na(countmat)] <- 0.5
    for(i in 1:nrow(countmat)) countmat[i,] <- 1000000*countmat[i,]/sum(countmat[i,])
    countmat <- log2(countmat)
    COUNTMAT <- countmat
    file <- paste0(dirout,"data_",welltype,"_",cutoff,".RData")
    save(CTRLMAT,PLATEMAT,COUNTMAT,file=file)
  }
}
