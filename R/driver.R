library(openxlsx)
library(tidyr)
library(stringr)
library(stringi)
library(reshape2)
#--------------------------------------------------------------------------------------
#' Code to run all signature concentration-response calculations
#'
#' @param dataset Name of the data set, produced by buildFCMAT2
#' @param sigcatalog Name of the signature catalog
#' @param sigset Name if the signature set. THis corresponds to a column in the signature catalog file
#' @param nullset The name of the NULL set if it is custom, default is NULL
#' @param nrandom.chems Number of random chemicals for the NULL distribution calculation, default is 1000
#' @param normfactor Normalization factor for the conc-reap plots, default is 7500
#' @param mc.cores Number of cores for paralell processing. Only works under Linux
#' @param bmr_scale Scaling factor from the NULL SD to BMD, default is 1.349
#' @param pval Threshold for cutoff distribution confidence interval. Default=0.05 indicates a 95% CI for the baseline distribution
#' @param nlowconc Only include the lowest nlowconc concentrations for each chemical
#' @param hccut The threshold for sigantures to be called a hit, default=0.95,
#' @param tccut The threshold for top/cutoff o be a hit, default =1.5
#' @param plotrange The concentration range for the conc-resp plots in uM, default is c(0.0001,100),
#' @param method signature scoring method in c("fc", "gsva", "gsea"), default is fc
#' @param celltype Name of the cull type, e.g. MCF7
#' @param do.build.random If TRUE, build the random dataset
#' @param do.run.random If True, run the calculations on the random data set. This is used
#'  to generate the NULL distribution
#' @param do.run.all If true, run the calculations on the real data set
#' @param do.scr.plots If TRUE, generate the signature concentration response plots
#' @param do.signature.summary.plot if TRUE, generate the summary plots
#' @param do.signature.pod If TRUE, generate the signature PODs
#' @param do.signature.pod.laneplot If TRUE, generate the signature lane plots (only useful for small sets of chemicals)
#' @param do.supertarget.boxplot If TRUE, generate the super target box plots
#' @param do.all If TRUE, do all steps from do.build.random to the end
#'
#'
#' Available data sets
#' * heparg2d_toxcast_pfas_pe1_normal
#' * mcf7_ph1_pe1_normal_block_123_allPG
#' * mcf7_ph1_pe1_normal_block_123_excludePG
#' * u2os_toxcast_pfas_pe1_normal
#' * PFAS_HepaRG
#' * PFAS_U2OS
#' * u2os_pilot_pe1_normal_null_pilot_lowconc
#' * u2os_toxcast_pfas_pe1_normal_refchems
#' * heparg2d_toxcast_pfas_pe1_normal_refchems
#' * DMEM_6hr_pilot_normal_pe_1 - MCF7 pilot
#' @export
#--------------------------------------------------------------------------------------
driver <- function(dataset="heparg2d_toxcast_pfas_pe1_normal_refchems",
                   sigcatalog="signatureDB_master_catalog 2021-05-10",
                   sigset="screen_large",
                   nullset=NULL,
                   nrandom.chems=1000,
                   normfactor=7500,
                   mc.cores=25,
                   bmr_scale=1.349,
                   pval=0.05,
                   nlowconc=2,
                   hccut=0.95,
                   tccut=1.5,
                   plotrange=c(0.0001,100),
                   method="fc",
                   celltype="HepaRG",
                   # do.build.random=T,
                   # do.run.random=T,
                   do.calc.cutoff=F,
                   do.run.all=T,
                   do.scr.plots=T,
                   do.signature.pod=T,
                   do.supertarget.boxplot=T,
                   do.all=F) {
  printCurrentFunction(paste(dataset,":",sigset))

  if(do.calc.cutoff) {
    cutoffCalc(dataset=dataset,
               sigcatalog=sigcatalog,
               sigset=sigset,
               method=method,
               pval=pval,
               nlowconc=nlowconc,
               mc.cores=mc.cores,
               do.load=T,do.cov=T,do.compare=F,to.file=F,verbose=F)
  }
  # if(do.build.random) {
  #   randomdata(dataset=dataset, nchem=nrandom.chems)
  # }
  #
  if(is.null(nullset)) nullset = paste0(dataset,"_RAND",nrandom.chems)
  # if(do.run.random || do.all){
  #   runAllSignatureCR(dataset=nullset,
  #                     nullset=nullset,
  #                     sigset=sigset,
  #                     sigcatalog=sigcatalog,
  #                     method = method,
  #                     normfactor=normfactor,
  #                     bmr_scale=bmr_scale,
  #                     do.plot = F,
  #                     do.cr=F,
  #                     mc.cores = c(mc.cores,mc.cores))
  # }
  if(do.run.all || do.all){
    runAllSignatureCR(dataset=dataset,
                      nullset=nullset,
                      sigset=sigset,
                      sigcatalog=sigcatalog,
                      method = method,
                      bmr_scale=bmr_scale,
                      normfactor=normfactor,
                      pval=pval,
                      nlowconc=nlowconc,
                      do.plot=T,
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
                                 pval=pval,
                                 plotrange=plotrange)
  }
  if(do.signature.pod || do.all) {
    signaturePOD(do.load=T,
                 sigset=sigset,
                 dataset=dataset,
                 method=method,
                 bmr_scale=bmr_scale,
                 hccut=hccut)
  }
  if(do.supertarget.boxplot) {
    superTargetBoxplot(to.file=T,
                       do.load=F,
                       dataset=dataset,
                       sigset=sigset,
                       method=method,
                       celltype=celltype,
                       hccut=hccut,
                       tccut=tccut,
                       cutoff=5)

    superTargetPODplot(to.file=T,
                       dataset=dataset,
                       sigset=sigset,
                       method=method,
                       celltype=celltype,
                       hccut=hccut,
                       tccut=tccut,
                       cutoff=5)

    superTargetStats(do.load=T,
                    dataset=dataset,
                    sigset=sigset,
                    method=method,
                    celltype=celltype,
                    hccut=hccut,
                    tccut=tccut,
                    cutoff=5)
  }
}
