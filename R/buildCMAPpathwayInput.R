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

  n <- nrow(mat)

  mat$index <- seq(from=1,to=n)

  name.list <- c("pathway","pathset","pathway_class","path_info" ,"gene_list","ngene")
  res.up <- as.data.frame(matrix(nrow=n,ncol=length(name.list)))
  names(res.up) <- name.list
  res.dn <- res.up

  pnames.up <- paste("CMAP_up",mat$trt_id,mat$index,sep="_")
  pnames.dn <- paste("CMAP_dn",mat$trt_id,mat$index,sep="_")
  res.up$pathway <- pnames.up
  res.dn$pathway <- pnames.dn

  res.up$pathset <- "CMAP"
  res.dn$pathset <- "CMAP"

  res.up$pathway_class <- NA
  res.dn$pathway_class <- NA

  pinfo <- paste(mat$name,mat$dsstox_sid,sep=" ")
  res.up$path_info <- pinfo
  res.dn$path_info <- pinfo

  res.up$gene_list <- mat$genes_up
  res.dn$gene_list <- mat$genes_dn

  res.up$ngene <- mat$sig_size
  res.dn$ngene <- mat$sig_size

  res <- rbind(res.up,res.dn)
  print(dim(res))
  res <- unique(res)
  print(dim(res))
  browser()
  CMAP_PATHWAYS <- res
  file <- "../input/processed_pathway_data/CMAP_PATHWAYS.RData"
  save(CMAP_PATHWAYS,file=file)
}
