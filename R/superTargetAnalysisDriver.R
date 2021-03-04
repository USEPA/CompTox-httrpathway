#--------------------------------------------------------------------------------------
#' Run analyses at the super_target level
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#' u2os_pilot_pe1_normal_null_pilot_lowconc
#' u2os_toxcast_pfas_pe1_normal_refchems
#' heparg2d_toxcast_pfas_pe1_normal_refchems
#'
#' After running this, run the following ...
#' superTargetPODplot
#' superTargetStats
#--------------------------------------------------------------------------------------
superTargetAnalysisDriver <- function(dataset="PFAS_U2OS",
                                      sigset="screen_large",
                                      method="fc",
                                      celltype="U2OS",
                                      hccut=0.95,
                                      tccut=1.5,
                                      minconc=0.0001,
                                      maxconc=100) {
  printCurrentFunction(paste(dataset,sigset,method))

  # input: SIGNATURE_CR_...RData
  # output: super_target_summary_ ...xlsx
  # output: super_target_boxplot_"..._all.RData"
  superTargetBoxplot(to.file=T,do.load=T,
                     dataset,sigset,method,celltype,
                     hccut,tccut,
                     minconc,maxconc)


  # input: SIGNATURE_CR_...RData
  # input super_target_boxplot_....xlsx"
  # output: super_target_summary_ ...xlsx

  superTargetSummary(do.load=F,verbose=T,
                     dataset,sigset,method,celltype,
                     hccut,tccut)

  #input: super_target_boxplot_..._all.RData"
  # output:super_target_boxplot_..._stats.xlsx"
  superTargetStats(do.load=T,
                   dataset,sigset,method,celltype,
                   hccut,tccut)

  # input: super_target_boxplot_...._summary.xlsx"
  # output: plots
  superTargetPODplot(to.file=T,
                     dataset,sigset,method,celltype,
                     hccut,tccut)

  # input: SIGNATURE_CR_...RData
  # output: plots
  superTargetOnOffTarget(to.file=T,do.load=F,
                         dataset,sigset,method,celltype,
                         hccut,tccut,
                         minconc,maxconc)

}

