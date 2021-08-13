#--------------------------------------------------------------------------------------
#' Merge multiple data sets for the same chemicals but across sigsets
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#' u2os_pilot_pe1_normal_null_pilot_lowconc
#' u2os_toxcast_pfas_pe1_normal_refchems
#' heparg2d_toxcast_pfas_pe1_normal_refchems
#--------------------------------------------------------------------------------------
signatureConcRespSummaryMerge <- function(do.load=T,
                                          dataset="mcf7_ph1_pe1_normal_block_123",
                                          sigset="screen_large",
                                          ssnew = c("stress_test","dorothea"),
                                          method="fc") {
  printCurrentFunction(paste(dataset,sigset,method))
  file0 = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
  print(file0)
  if(do.load) {
    load(file=file0)
    mat = SIGNATURE_CR
    MAT <<- mat
  }
  mat = MAT
  if(!is.element("bmr_scale",names(mat))) {
    mat$bmr_scale = 1.349
  }
  cat(nrow(mat),"\n")
  for(ss in ssnew) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",ss,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    temp = SIGNATURE_CR
    mat = rbind(mat,temp)
    cat(nrow(mat),"\n")
  }

  SIGNATURE_CR = mat
  save(SIGNATURE_CR,file=file0)
}

