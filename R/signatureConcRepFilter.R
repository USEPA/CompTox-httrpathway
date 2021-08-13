library(openxlsx)
library(tidyr)
library(stringr)
library(stringi)
library(reshape2)
library(gplots)
library(RColorBrewer)
#--------------------------------------------------------------------------------------
#' Filter the conc-repons data
#' @param method signature scoring method in c("fc", "gsva", "mygsea")
#'
#'
#' Error bars are exp(er)*qt(.025,4) = exp(er)*2.7765
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' mcf7_ph1_pe1_normal_block_123_allPG
#' u2os_toxcast_pfas_pe1_normal
#--------------------------------------------------------------------------------------
signatureConcRepFilter <- function(to.file=F,
                                   do.plot=F,
                                   do.load=T,
                                   hccut=0.9,
                                   tccut=1.5,
                                   dataset="heparg2d_toxcast_pfas_pe1_normal",
                                   sigset="screen_large",
                                   method="fc",
                                   do.pfas=F) {
  printCurrentFunction(paste(dataset,sigset,method))

  if(do.load) {
    cat("do.load\n")
    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat <- SIGNATURE_CR
    MAT <<- mat
  }

  #######################################################################################
  # Do all chemicals
  #######################################################################################
  mat <- MAT
  cat(nrow(mat),"\n")
  mat <- mat[mat$top_over_cutoff>tccut,]
  cat(nrow(mat),"\n")
  mat <- mat[mat$hitcall>hccut,]
  cat(nrow(mat),"\n")
  mat <- mat[order(mat$top_over_cutoff,decreasing=T),]
  cat(nrow(mat),"\n")
  mat$proper_name <- mat$name
  file <- paste0("../output/signature_conc_resp_filtered/signature_conc_resp_filtered ",dataset,"_",sigset,"_",hccut,"_",tccut,".RData")
  SIGNATURE_CR <- mat
  save(SIGNATURE_CR,file=file)
  if(do.plot) {
    if(to.file) {
      fname <- paste0("../output/signature_conc_resp_filtered/signature_conc_resp_filtered ",dataset,"_",sigset,"_",hccut,"_",tccut,".pdf")
      pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    }
    par(mfrow=c(3,2),mar=c(4,4,2,2))
    for(i in 1:nrow(mat)){
      if(i%%1000==0) cat(i," out of ",nrow(mat),"\n")
      tryCatch({
        signatureConcRespPlot(mat[i,])
        if(!to.file) browser()
      }, warning = function(w) {
        cat("WARNING:\n")
      }, error = function(e) {
        cat("ERROR\n")
      })
    }

    if(!to.file) browser()
    else dev.off()
  }

  if(do.pfas) {
    #######################################################################################
    # Do PFAS chemicals
    #######################################################################################
    file <- "../input/pfas/QC sample map 2020-05-04.xlsx"
    qc <- read.xlsx(file)
    #qc <- qc[!is.element(qc$score,c("H","M")),]
    spid.list <- qc$spid
    mat <- mat[is.element(mat$sample_id,spid.list),]
    if(nrow(mat)>0) {
      if(to.file) {
        fname <- paste0("../output/signature_conc_resp_filtered/signature_conc_resp_filtered pfas ",dataset,"_",sigset,"_",hccut,"_",tccut,".pdf")
        pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
      }
      par(mfrow=c(3,2),mar=c(4,4,2,2))
      file <- paste0("../output/signature_conc_resp_fsignatureConcRepFilteriltered/signature_conc_resp_filtered pfas ",dataset,"_",sigset,"_",hccut,"_",tccut,".RData")
      SIGNATURE_CR <- mat
      save(SIGNATURE_CR,file=file)

      for(i in 1:nrow(mat)){
        cat(i," out of ",nrow(mat),"\n")
        tryCatch({
          signatureConcRespPlot(mat[i,])
          if(!to.file) browser()
        }, warning = function(w) {
          cat("WARNING:\n")
        }, error = function(e) {
          cat("ERROR\n")
        })
      }
      signatureConcRepFilter
      if(!to.file) browser()
      else dev.off()
    }
  }
}
