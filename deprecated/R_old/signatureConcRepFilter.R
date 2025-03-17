#--------------------------------------------------------------------------------------
#' Filter the conc-repons data
#' @param to.file If TRUE, send plots to a file
#' @param do.plot If TRUE do the plotting
#' @param do.load If TRUE, load the data file
#' @param hccut Exclude rows with hitcall below this value
#' @param tccut Exclude rows with top_over_cutoff below this value
#' @param dataset Dataset to use
#' @param sigset Signature set to use
#' @param method signature scoring method in c("fc", "gsva", "gsea")
#' @param do.pfas=F, if T handle pfas chemicals
#'
#' Error bars are exp(er)*qt(.025,4) = exp(er)*2.7765
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' mcf7_ph1_pe1_normal_block_123_allPG
#' u2os_toxcast_pfas_pe1_normal
#' @importFrom grDevices pdf dev.off
#' @importFrom openxlsx read.xlsx
#' @importFrom graphics par
#' @export signatureConcRepFilter
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
    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RDS")
    print(file)
    MAT <<- readRDS(file)
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
  file <- paste0("../output/signature_conc_resp_filtered/signature_conc_resp_filtered ",dataset,"_",sigset,"_",hccut,"_",tccut,".RDS")
  SIGNATURE_CR <- mat
  saveRDS(SIGNATURE_CR,file)
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
      file <- paste0("../output/signature_conc_resp_fsignatureConcRepFilteriltered/signature_conc_resp_filtered pfas ",dataset,"_",sigset,"_",hccut,"_",tccut,".RDS")
      
      SIGNATURE_CR <- mat
      saveRDS(SIGNATURE_CR,file)

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
      signatureConcRepFilter #??
      if(!to.file) browser()
      else dev.off()
    }
  }
}
