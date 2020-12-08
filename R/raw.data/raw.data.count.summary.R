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
raw.data.count.summary <- function(to.file=F,
                                   dirin="../../../Raw/ctrl/",
                                   dirout="analysis/raw_count_analysis/data/") {
  printCurrentFunction()
  if(to.file) {
    fname <- paste0(dirout,"raw_data_count_summary.pdf")
    pdf(file=fname,width=7,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,2,2))

  file.list <- list.files(dirin)
  for(f in file.list) {
    file <- paste0(dirin,f)
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
    countmat[is.na(countmat)] <- 0.5
    countmat[countmat<10] <- 0
    rs <- log10(rowSums(countmat))
    hist(rs,main=f,xlab="log10(total well counts)")
    if(!to.file) browser()
  }
  if(to.file) dev.off()
}
