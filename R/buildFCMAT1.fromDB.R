library(tidyverse)
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
#' @param pg.filter.file An optional file to use in filtering out bad plate groups
#' @param do.load If TRUE, read the large input data file into memory
#' @return A file with the FCMAT1 data is written to "../input/fcdata/FCMAT1_",dataset,".RData"
#'
#' @export
#--------------------------------------------------------------------------------------
buildFCMAT1.fromDB <- function(dataset="tox21_cpp5_u2os_pe1_normal",
                               dir="../input/fcdata/new_versions/",
                               infile="httr_tox21_cpp5_u2os_FCmat1-meanncnt0_5-plateteffect_1-shrinkage_normal.RData",
                               pg.filter.file=NULL,
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

  temp = unique(FCMAT1[,c("chem_id","chem_name","dtxsid","pg_id","block_id")])
  names(temp) = c("sample_id","name","dtxsid","pg_id","block_id")
  file = paste0(dir,"PG_MAP_",dataset,".xlsx")
  write.xlsx(temp,file=file)
  browser()

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
  file <- paste0("../input/fcdata/FCMAT1_",dataset,".RData")
  cat("start saving\n")
  save(FCMAT1,file=file)
  cat("finish saving RData file\n")
}
