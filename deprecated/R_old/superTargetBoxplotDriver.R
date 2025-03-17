#--------------------------------------------------------------------------------------
#' Do  series of runs of superTargetBoxplotDriver
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#' mcf7_ph1_pe1_normal_all_pg
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#'
#--------------------------------------------------------------------------------------
superTargetBoxplotDriver <- function(dataset="mcf7_ph1_pe1_normal_block_123",
                                     sigset="screen_large",
                                     method="fc",
                                     celltype="MCF7",
                                     tcset=c(1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0),
                                     hcset=c(0.6,0.8,0.9)){
  printCurrentFunction(paste(dataset,sigset,method))

  do.load=T
  for(tccut in tcset) {
    for(hccut in hcset) {
      superTargetBoxplot(to.file=T,
                         do.load=do.load,
                         dataset=dataset,
                         sigset=sigset,
                         method=method,
                         celltype=celltype,
                         hccut=hccut,
                         tccut=tccut)
      do.load=F
    }
  }


}

