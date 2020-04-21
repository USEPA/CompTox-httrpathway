#--------------------------------------------------------------------------------------
#' Build the standard input file for the Bioplanet signatures
#'
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
signatureBuildBioplanet <- function(){
  printCurrentFunction()

  file <- "../input/signatures/bioplanet_PATHWAYS.RData"
  load(file)

  mat <- bioplanet_PATHWAYS

  names(mat) <- c("signature","source","class","description","gene.list","ngene")
  mat$parent <- mat$signature
  mat$source <- "Bioplanet"
  mat$type <- "nondirectional"
  mat$direction <- "nondirectional"
  mat$subsource <- "-"
  name.list <- c("signature","parent","source","type","direction","description","subsource","ngene","gene.list")
  mat <- mat[,name.list]

  Bioplanet_signatures <- mat
  file <- "../input/signatures/Bioplanet_signatures.RData"
  save(Bioplanet_signatures,file=file)
}
