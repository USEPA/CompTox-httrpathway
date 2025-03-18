#--------------------------------------------------------------------------------------
#'
#' Export the signature-wise cutoffs
#' @param signature_cr dataframe returned by the signatureConcResp function
#' @export exportSignatureCutoffs
#--------------------------------------------------------------------------------------
exportSignatureCutoffs <- function(signature_cr) {
  printCurrentFunction()
  signature_cutoff = unique(signature_cr[,c("signature","signature_size","cutoff")])
  return(signature_cutoff)
}
