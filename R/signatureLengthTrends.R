library(openxlsx)
library(tidyr)
library(stringr)
library(stringi)
library(reshape2)
library(gplots)
library(RColorBrewer)
#--------------------------------------------------------------------------------------
#' Calulate trends in conc-response data as function of signature length
#' @param method signature scoring method in c("fc", "gsva", "mygsea")
#'
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#'  u2os_toxcast_pfas_pe1_normal
#'
#' Error bars are exp(er)*qt(.025,4) = exp(er)*2.7765
#--------------------------------------------------------------------------------------
signatureLengthTrends <- function(to.file=F,
                                  do.load=F,
                                  dataset="u2os_toxcast_pfas_pe1_normal",
                                  sigcatalog="signatureDB_master_catalog 2020-08-14",
                                  sigset="screen_large",
                                  method="fc",
                                  celltype="U2OS") {
  printCurrentFunction(paste(dataset,sigset,method))
  if(to.file) {
    fname <- paste0("../output/signature_length/signature_length_",celltype,
                    "_",dataset,"_",sigset,"_",method,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,2,2))

  if(do.load) {
    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat <- SIGNATURE_CR
    MAT <<- mat
  }
  mat <- MAT
  mat <- mat[mat$hitcall>0.9,]
  mat$auc <- abs(mat$top) *(3-log10(mat$bmd))
  mat$ci <- exp(mat$er)*2.7765
  vals <- c(0,20,40,60,80,100,200,300,400,500,1000)
  nval <- length(vals)
  slen <- NULL
  tc <- NULL
  hc <- NULL
  bmd <- NULL
  auc <- NULL
  ci <- NULL
  for(i in 1:nval) {
    bottom <- vals[i]
    if(i<nval) top <- vals[i+1]
    else top=1100

    temp <- mat[mat$signature_size>=bottom,]
    temp <- temp[temp$signature_size<top,]
    n <- nrow(temp)
    xx <- vector(mode="character",length=n)
    xx[] <- paste0(bottom,"-",top)
    slen <- c(slen,xx)
    yy <- temp$hitcall
    hc <- c(hc,yy)
    yy <- temp$top_over_cutoff
    tc <- c(tc,yy)
    yy <- temp$bmd
    bmd <- c(bmd,yy)
    yy <- temp$auc
    auc <- c(auc,yy)
    yy <- temp$ci
    ci <- c(ci,yy)
  }
  boxplot(hc~slen,cex.lab=1.2,cex.axis=1.2,xlab="signature length",ylab="Hitcall",main=paste("Hitcall",celltype))
  boxplot(tc~slen,cex.lab=1.2,cex.axis=1.2,ylim=c(0,10),xlab="signature length",ylab="T/C",main=paste("T/C",celltype))
  boxplot(bmd~slen,log="y",cex.lab=1.2,cex.axis=1.2,xlab="signature length",ylab="BMD",main=paste("BMD",celltype))
  boxplot(auc~slen,cex.lab=1.2,cex.axis=1.2,ylim=c(0,2),xlab="signature length",ylab="AUC",main=paste("AUC",celltype))
  boxplot(ci~slen,cex.lab=1.2,cex.axis=1.2,ylim=c(0,0.2),xlab="signature length",ylab="CI",main=paste("CI",celltype))

  if(to.file) dev.off()
  else browser()
}
