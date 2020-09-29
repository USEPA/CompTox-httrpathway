library(openxlsx)
library(tidyr)
library(stringr)
library(stringi)
library(reshape2)
#--------------------------------------------------------------------------------------
#' Code to run all calculations
#' @param method signature scoring method in c("fc", "gsva", "mygsea")
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#' mcf7_ph1_pe1_normal_all_pg
#' u2os_toxcast_pfas_pe1_normal
#' DMEM_6hr_pilot_normal_pe_1
#' mcf7_ph1_pe1_normal_all_pg
#' mcf7_ph1_pe1_normal_block_123
#'
#'  signatureDB_wgcna_mcf7_ph1_pe1_normal_good_pg_MCF7_12_10_catalog
#'  signatureDB_master_catalog 2020-07-10 / screen_large
#'
#'
#'
#--------------------------------------------------------------------------------------
driver <- function(dataset="mcf7_ph1_pe1_normal_block_123",
                   sigcatalog="signatureDB_master_catalog 2020-09-16",
                   sigset="screen_large",
                   nullset=NULL,
                   nrandom.chems=1000,
                   normfactor=7500,
                   mc.cores=20,
                   method="fc",
                   do.build.random=T,
                   do.run.random=T,
                   do.run.all=T,
                   do.scr.plots=T,
                   do.signature.summary.plot=T,
                   do.signature.pod=T,
                   do.signature.pod.laneplot=F,
                   do.all=F) {
  printCurrentFunction(paste(dataset,":",sigset))

  if(do.build.random) {
    randomdata(dataset=dataset, nchem=nrandom.chems)
  }

  if(is.null(nullset)) nullset <- paste0(dataset,"_RAND",nrandom.chems)
  if(do.run.random || do.all){
    runAllSignatureCR(dataset=nullset,
                      nullset=nullset,
                      sigset=sigset,
                      sigcatalog=sigcatalog,
                      method = method,
                      normfactor=normfactor,
                      do.plot = F,
                      do.cr=F,
                      mc.cores = c(mc.cores,mc.cores))
  }
  if(do.run.all || do.all){
    runAllSignatureCR(dataset=dataset,
                      nullset=nullset,
                      sigset=sigset,
                      sigcatalog=sigcatalog,
                      method = method,
                      normfactor=normfactor,
                      do.plot = T,
                      do.cr=T,
                      mc.cores = c(mc.cores,mc.cores))
    cat("Look for output in \n
        ../output/signature_score_summary/\n
        ../output/signature_conc_resp_plots/\n
        ../output/signature_conc_resp_summary/\n
        \n")
  }
  if(do.all || do.scr.plots){
    signatureConcRespPlotWrapper(sigset=sigset,
                                  dataset=dataset,
                                  sigcatalog=sigcatalog,
                                  method=method,
                                  mc.cores=mc.cores,
                                  do.load=T,
                                  pval = .05,
                                  nametag = NULL)
  }
  if(do.signature.summary.plot || do.all) {
    signatureClassSummaryPlot(to.file=T,dataset=dataset,
                              sigcatalog=sigcatalog,
                              sigset=sigset,
                              method = method)
    #signatureClassSummaryDotPlot(to.file=T,dataset=dataset,
    #                             sigcatalog=sigcatalog,
    #                            sigset=sigset,
    #                             method = method)
  }
  if(do.signature.pod || do.all) {
    signaturePOD(do.load=T,
                 sigset=sigset,
                 dataset=dataset,
                 method=method,
                 hccut=0.9)
  }
  if(do.signature.pod.laneplot) {
    podLaneplot(to.file=T,
                dataset=dataset,
                sigset=sigset,
                method=method,
                hccut=0.9,
                plot.signature_min=F)
  }
}
