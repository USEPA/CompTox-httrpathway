#--------------------------------------------------------------------------------------
#'
#' Export the signature-wise cutoffs
#' @param do.load If TRUE, load hte large data file
#' @param dataset The name of the HTTr data set to use
#' @param sigset The name of the signature set to use
#' @param THe scoring method to use
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123_allPG
#' mcf7_ph1_pe1_normal_block_123_excludePG
#' u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#' u2os_pilot_pe1_normal_null_pilot_lowconc
#'
#--------------------------------------------------------------------------------------
exportSignatureCutoffs <- function(do.load=F,
                                   dataset="mcf7_ph1_pe1_normal_block_123_excludePG",
                                   sigset="screen_large",
                                   method="fc") {
  printCurrentFunction()
  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    MAT <<- mat
  }
  mat = MAT
  temp = unique(mat[,c("signature","signature_size","cutoff")])
  file = paste0("../output/signature_cutoff/signature_cutoff_",sigset,"_",dataset,"_",method,"_0.05.xlsx")
  write.xlsx(temp,file)
}
