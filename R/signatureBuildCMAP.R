#--------------------------------------------------------------------------------------
#' Build the standard input file for the CMAP signatures
#'
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
signatureBuildCMAP <- function(){
  printCurrentFunction()

  file <- "../input/signatures/cmap-sig-2019-11-06-a.tsv"
  mat <- read.table(file,header=T,sep="\t",quote="",stringsAsFactors=F)
  x <- paste("CMAP",mat[,3],mat[,5],mat[,6])
  mat[,1] <- x
  x <- mat[,1]
  y <- seq(from=1,to=length(x))
  x <- paste(x,y)
  mat[,1] <- x
  mat.up <- mat[,c(1:6,7)]
  mat.dn <- mat[,c(1:6,8)]
  name.list <- c("signature","cell","name","dtxsid","conc","ngene","gene.list")
  names(mat.up) <- name.list
  names(mat.dn) <- name.list

  mat.up$signature <- paste0(mat.up$signature," ",mat.up$ngene)
  mat.dn$signature <- paste0(mat.dn$signature," ",mat.dn$ngene)
  mat.up$subsource <- "-"
  mat.dn$subsource <- "-"

  mat.up$parent <- mat.up$signature
  mat.dn$parent <- mat.dn$signature
  mat.up$type <- "directional"
  mat.dn$type <- "directional"
  mat.up$direction <- "up"
  mat.dn$direction <- "dn"
  mat.up$source <- "CMAP"
  mat.dn$source <- "CMAP"

  mat.up$signature <- paste0(mat.up$signature," up")
  mat.dn$signature <- paste0(mat.dn$signature," dn")

  mat.up$description <- paste(mat.up$name,mat.up$dtxsid)
  mat.dn$description <- paste(mat.dn$name,mat.dn$dtxsid)


  name.list <- c("signature","parent","source","type","direction","description","subsource","ngene","gene.list")
  mat.up <- mat.up[,name.list]
  mat.dn <- mat.dn[,name.list]
  mat<- rbind(mat.up,mat.dn)

  mat <- mat[order(mat$signature),]
  CMAP_signatures <- mat
  file <- "../input/signatures/CMAP_signatures.RData"
  save(CMAP_signatures,file=file)
}
