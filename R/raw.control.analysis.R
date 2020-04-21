#source("https://bioconductor.org/biocLite.R")
#.libPaths( "C:/Program Files/R/R-3.3.3/library")
library(fastmatch)
library(reshape2)
#source("R/general_utils.R")
#options(java.parameters = "-Xmx8000m")
#library(httRws)
library(openxlsx)
#--------------------------------------------------------------------------------------
#' Prepare the raw data for the replicated chemicals
#'
#' @param dir The directory to work in
#'
#--------------------------------------------------------------------------------------
raw.control.analysis <- function(to.file=F,
                                 do.read=F,
                                 do.prep=F,
                                 dirin="../../../Raw/ctrl/",
                                 welltype="dmso",
                                 dirout="analysis/raw_count_analysis/",
                                 cutoff=5) {
  printCurrentFunction()
  if(do.read) {
    file <- paste0(dirin,"httr-ph-i-raw-pl-bl_",welltype,"-v1.tsv")
    print(file)
    mat <- read.delim(file,stringsAsFactors=F)
    CTRLMAT <<- mat
    pmat <- unique(CTRLMAT[,c("plate_id","block_id","well_id","mapd_frac","n_reads")])
    PLATE.MAT <<- pmat
  }
  if(do.prep) {
    mat <- CTRLMAT
    pmat <- PLATE.MAT

    plate.list <- sort(unique(mat[,"plate_id"]))
    plate.list <- plate.list[!is.element(plate.list,"plate_id")]
    nplate <- length(plate.list)
    probe.list <- sort(unique(mat[,"probe_id"]))
    nprobe <- length(probe.list)
    res.A <- as.data.frame(matrix(nrow=nprobe,ncol=nplate))
    rownames(res.A) <- probe.list
    res.A[] <- 0
    res.B <- res.A
    for(i in 1:nplate) {
      plate <- plate.list[i]
      temp1 <- mat[is.element(mat[,"plate_id"],plate),]
      block_id <- unique(temp1$block_id)
      cat(plate,nrow(temp1),block_id,"\n")
      print(unique(temp1[,"well_id"]))
      colname.A <- paste0(plate,"_",block_id,"_A")
      colname.B <- paste0(plate,"_",block_id,"_B")
      names(res.A)[i] <- colname.A
      names(res.B)[i] <- colname.B
      x1 <- temp1[is.element(temp1[,"well_id"],"E01"),c("probe_id","probe_count")]
      x2 <- temp1[is.element(temp1[,"well_id"],"F01"),c("probe_id","probe_count")]
      x1 <- x1[order(x1$probe_id),]
      x2 <- x2[order(x2$probe_id),]
      res.A[x1$probe_id,colname.A] <- as.numeric(x1$probe_count)
      res.B[x2$probe_id,colname.B] <- as.numeric(x2$probe_count)
    }
    RES.A <<- res.A
    RES.B <<- res.B
  }
  res.A <- RES.A
  res.B <- RES.B
  if(to.file) {
    fname <- paste0(dirout,"raw_count_comparison_",welltype,"_residuals.pdf")
    pdf(file=fname,width=10,height=6,pointsize=12,bg="white",paper="USr",pagecentre=T)
  }
  par(mfrow=c(4,1),mar=c(2,2,4,4))
  res <- cbind(res.A,res.B)
  col.list <- names(res)
  x <- str_split(col.list,"_")
  colmap <- as.data.frame(matrix(unlist(x),nrow=length(col.list),byrow=T))
  colmap <- colmap[order(colmap[,2]),]
  col.list <- paste0(colmap[,1],"_",colmap[,2],"_",colmap[,3])

  ncol <- length(col.list)
  #ncol <- 20
  ncoli <- ncol
  #ncoli <- 100
  for(i in 1:ncoli) {
    labels <- NULL
    values <- NULL
    blocks <- NULL
    cat(i,"\n")
    for(j in i:ncol) {
      colA <- col.list[i]
      colB <- col.list[j]
      blockB <- str_split(colB,"_")[[1]][2]
      xA <- res[,colA]
      xB <- res[,colB]
      mask <- xA
      mask[] <- 1
      mask[xA<cutoff] <- 0
      mask[xB<cutoff] <- 0
      xA <- xA[mask==1]
      xB <- xB[mask==1]
      xA <- log2(xA)
      xB <- log2(xB)
      vtemp <- xA-xB
      ltemp <- vtemp
      btemp <- vtemp
      ltemp[] <- colB
      btemp[] <- as.numeric(blockB)
      labels <- c(labels,ltemp)
      values <- c(values,vtemp)
      blocks <- c(blocks,btemp)
      #fit <- lm(xB~xA+0)
      #r2 <- summary(fit)$adj.r.squared
      #doit <- F
      ##if(r2<0.95) doit <- T
      #if(i+j < 5) doit <- T
      #if(doit) {
      #  plot(xB~xA,pch=".",xlab=colA,ylab=colB,xlim=c(2,15),ylim=c(2,15))
      #  text(2,14,colA,pos=4,cex=0.9)
      #  text(2,12.5,colB,pos=4,cex=0.9)
      #  text(2,11,format(r2,digits=3),pos=4,cex=0.9)
      #  lines(c(0,100),c(0,100))
      #  if(!to.file) browser()
      #}
    }
    boxplot(values~labels,main=colA,ylim=c(-6,6),outline=F)
    lines(c(-1000,1000),c(0,0),lwd=2)
    if(!to.file) browser()
  }
  if(to.file) dev.off()
}
