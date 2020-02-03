#--------------------------------------------------------------------------------------
#' Build the standard input file for the CMAP signatures
#'
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
buildCMAPpathwayInput = function(){
  printCurrentFunction()

  file <- "../input/processed_pathway_data/cmap-sig-2019-11-06-a.tsv"
  mat <- read.table(file,header=T,sep="\t",quote="",stringsAsFactors=F)
  mat.up <- mat[,c(1:6,7)]
  mat.dn <- mat[,c(1:6,8)]
  name.list <- c("signature","cell","name","dtxsid","conc","ngene","gene.list")
  names(mat.up) <- name.list
  names(mat.dn) <- name.list

  mat.up$signature.parent <- mat.up$signature
  mat.dn$signature.parent <- mat.dn$signature
  mat.up$signature.type <- "bidirectional"
  mat.dn$signature.type <- "bidirectional"
  mat.up$direction <- "up"
  mat.dn$direction <- "dn"
  mat.up$signature.set <- "CMAP"
  mat.dn$signature.set <- "CMAP"

  mat.up$signature <- paste0(mat.up$signature,"_up")
  mat.dn$signature <- paste0(mat.dn$signature,"_dn")

  mat.up$description <- paste(mat.up$name,mat.up$dtxsid,mat.up$conc,mat.up$cell)
  mat.dn$description <- paste(mat.dn$name,mat.dn$dtxsid,mat.dn$conc,mat.dn$cell)


  name.list <- c("signature","parent","set","type","direction","description","ngene","gene.list")
  mat.up <- mat.up[,name.list]
  mat.dn <- mat.dn[,name.list]
  mat<- rbind(mat.up,mat.dn)
  browser()

  genes <- paste0(mat$genes_up,"|",mat$genes_dn)
  mat <- cbind(mat,genes)
  mat <- mat[,c(1:6,9)]
  n <- nrow(mat)

  mat$index <- seq(from=1,to=n)

  name.list <- c("pathway","pathset","pathway_class","path_info" ,"gene_list","ngene")
  res <- as.data.frame(matrix(nrow=n,ncol=length(name.list)))
  names(res) <- name.list

  pnames <- paste("CMAP",mat$trt_id,mat$index,sep="_")
  res$pathway <- pnames

  res$pathset <- "CMAP"
  res$pathway_class <- NA

  pinfo <- paste(mat$name,mat$dsstox_sid,sep=" ")
  res$path_info <- pinfo

  res$gene_list <- mat$genes

  res$ngene <- mat$sig_size

  print(dim(res))
  res <- unique(res)
  print(dim(res))
  browser()
  CMAP_PATHWAYS <- res
  file <- "../input/processed_pathway_data/CMAP_PATHWAYS.RData"
  save(CMAP_PATHWAYS,file=file)
}
