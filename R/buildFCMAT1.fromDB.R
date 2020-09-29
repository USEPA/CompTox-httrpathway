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
#' httr_mcf7_ph1_FCmat1_meanncnt0_5-plateteffect_1-shrinkage_normal
library(tidyverse)
#--------------------------------------------------------------------------------------
buildFCMAT1.fromDB <- function(dataset="mcf7_ph1_pe1_normal_block_123",
                               dir="../input/fcdata/new_versions/",
                               infile="httr_mcf7_ph1_FCmat1_meanncnt0_5-plateteffect_1-shrinkage_normal_good_pg.RData",
                               pg.filter.file="httr_mcf7_ph1_flagged_pg_block_123.xlsx",
                               do.load=T){
  printCurrentFunction()
  if(do.load) {
    cat("   start loading RData file\n")
    infile <- paste0(dir,infile)
    print(infile)
    load(file=infile)
    cat("finished loading\n")
    FCMAT1.0 <<- FCMAT1
  }
  FCMAT1 = FCMAT1.0
  cat("initial:",dataset,":",nrow(FCMAT1),"\n")

  filter = NULL
  if(!is.null(pg.filter.file)) {
    file = paste0("../input/fcdata/",pg.filter.file)
    filter = read.xlsx(file)
    filter$pg_id = as.character(filter$pg_id)
  }
  if(!is.null(filter)) {
    good.pg = filter[is.element(filter$pg_flag,"OK"),"pg_id"]
    FCMAT1 = FCMAT1[is.element(FCMAT1$pg_id,good.pg),]
  }
  name.list <- names(FCMAT1)
  name.list[is.element(name.list,"gene_symbol")] <- "gene"
  names(FCMAT1) <- name.list
  cat("final:",dataset,":",nrow(FCMAT1),"\n")
  #browser()
  file <- paste0("../input/fcdata/FCMAT1_",dataset,".RData")
  cat("start saving\n")
  save(FCMAT1,file=file)
  cat("finish saving RData file\n")
}
