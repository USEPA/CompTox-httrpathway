library(openxlsx)
library(tidyr)
library(stringr)
library(stringi)
library(reshape2)
#--------------------------------------------------------------------------------------
#' Code to run all calculations
#' @param method signature scoring method in c("fc", "gsva", "mygsea")
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#' u2os_pilot_pe1_normal_null_pilot_lowconc
#'
#' DMEM_6hr_pilot_normal_pe_1 - MCF7 pilot
#'
#' u2os_toxcast_pfas_pe1_normal_refchems
#' heparg2d_toxcast_pfas_pe1_normal_refchems
#'
#--------------------------------------------------------------------------------------
driver <- function(dataset="DMEM_6hr_pilot_normal_pe_1",
                   sigcatalog="signatureDB_master_catalog 2021-03-05",
                   sigset="screen_large",
                   nullset=NULL,
                   nrandom.chems=1000,
                   normfactor=7500,
                   mc.cores=20,
                   bmr_scale=1.349,
                   plotrange=c(0.0001,100),
                   method="fc",
                   celltype="MCF7",
                   do.build.random=T,
                   do.run.random=T,
                   do.run.all=T,
                   do.scr.plots=T,
                   do.signature.summary.plot=T,
                   do.signature.pod=T,
                   do.signature.pod.laneplot=T,
                   do.supertarget.boxplot=T,
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
                      bmr_scale=bmr_scale,
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
                      bmr_scale=bmr_scale,
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
                                 bmr_scale=bmr_scale,
                                 mc.cores=mc.cores,
                                 do.load=T,
                                 pval = .05,
                                 plotrange=plotrange)
  }
  if(do.signature.summary.plot || do.all) {
    signatureClassSummaryPlot(to.file=T,dataset=dataset,
                              sigcatalog=sigcatalog,
                              sigset=sigset,
                              method = method,
                              bmr_scale=bmr_scale)
  }
  if(do.signature.pod || do.all) {
    signaturePOD(do.load=T,
                 sigset=sigset,
                 dataset=dataset,
                 method=method,
                 bmr_scale=bmr_scale,
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
  if(do.supertarget.boxplot) {
    superTargetBoxplot(to.file=T,
                       do.load=F,
                       dataset=dataset,
                       sigset=sigset,
                       method=method,
                       celltype=celltype,
                       hccut=0.95,
                       tccut=1.5)
  }
}
