#' get the list of genes used in a data set
#'
#' @param dataset name ofthe data set to use
#' heparg2d_toxcast_pfas_pe1_normal
#' u2os_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
expressedGeneList = function(dataset="mcf7_ph1_pe1_normal_block_123") {
  printCurrentFunction(dataset)
  dir = "../input/fcdata/"
  file = paste0(dir,"FCMAT2_",dataset,".RData")
  print(file)
  load(file=file)
  genelist = colnames(FCMAT2)
  cat(length(genelist),"\n")
  file = paste0(dir,"genelist_",dataset,".xlsx")
  write.xlsx(genelist,file)
}

