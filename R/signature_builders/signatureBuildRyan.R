#--------------------------------------------------------------------------------------
#' Build the standard input file for the Ryan signatures
#'
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
signatureBuildRyan <- function(){
  printCurrentFunction()

  file <- "../input/signatures/RYAN_PATHWAYS.RData"
  load(file)
  mat <- RYAN_PATHWAYS
  mat[1,"pathway"] <- "RYAN_ESTROGEN_RECEPTOR_ALPHA_up"
  mat[2,"pathway"] <- "RYAN_ESTROGEN_RECEPTOR_ALPHA_dn"
  names(mat) <- c("signature","set","class","description","gene.list","ngene")
  mat$parent <- " RYAN_ESTROGEN_RECEPTOR_ALPHA" #mat$signature
  mat$source <- "Ryan"
  mat$type <- "directional"
  mat[1,"direction"] <- "up"
  mat[2,"direction"] <- "dn"
  mat$subsource <- "-"
  name.list <- c("signature","parent","source","type","direction","description","subsource","ngene","gene.list")
  mat <- mat[,name.list]
  Ryan_signatures <- mat
  file <- "../input/signatures/Ryan_signatures.RData"
  save(Ryan_signatures,file=file)
}
