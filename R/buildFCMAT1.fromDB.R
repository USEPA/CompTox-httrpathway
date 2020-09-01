#--------------------------------------------------------------------------------------
#' Build the FCMAT1 data set
#'
#' version to start with Logan's database export
#' The difference between this version and the original is that there are extra columns
#' The function just changes one column name and writes the file to a standard name and place
#'
#' @param dataset The name to give to the data set
#' @param dir The directory from which to read all of the raw files
#' @param infile The nae of the input file
#' @return A file with the FCMAT1 data is written to "../input/fcdata/FCMAT1_",dataset,".RData"
#'
#' @export
library(tidyverse)
#--------------------------------------------------------------------------------------
buildFCMAT1.fromDB <- function(dataset="u2os_toxcast_pfas_pe1_normal",
                               dir="../input/fcdata/new_versions/",
                               infile="httr_u2os_toxcast_pfas_FCmat1-meanncnt0_5-plateteffect_1-shrinkage_normal.RData"){
  printCurrentFunction()
  cat("   start loading RData file\n")
  infile <- paste0(dir,infile)
  print(infile)
  load(file=infile)
  cat("finished loading\n")
  name.list <- names(FCMAT1)
  name.list[is.element(name.list,"gene_symbol")] <- "gene"
  names(FCMAT1) <- name.list
  file <- paste0("../input/fcdata/FCMAT1_",dataset,".RData")
  cat("start saving\n")
  save(FCMAT1,file=file)
  cat("   start saving RData file\n")
  save(FCMAT1,file=file)
  cat("finish saving RData file\n")
}
