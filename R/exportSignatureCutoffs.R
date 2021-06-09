#--------------------------------------------------------------------------------------
#'
#' Export the signature-wise cutoffs
#' @param do.load If TRUE, load hte large data file
#' @param dataset The name of the HTTr data set to use
#' @param sigset The name of the signature set to use
#' @param method The scoring method to use
#' @export
#--------------------------------------------------------------------------------------
exportSignatureCutoffs <- function(do.load=F,
                                   dataset="heparg2d_toxcast_pfas_pe1_normal",
                                   sigset="screen_large",
                                   method="fc") {
  printCurrentFunction()
  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/old/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
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
