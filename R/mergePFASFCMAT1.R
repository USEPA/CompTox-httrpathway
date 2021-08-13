library(tidyverse)
#--------------------------------------------------------------------------------------
#' special code the merge the PFAS replacemnte data in with the earlier data fro U2OS and HepaRG
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
#' * heparg2d_toxcast_pfas_pe1_normal
#' * u2os_toxcast_pfas_pe1_normal
#' @export
#--------------------------------------------------------------------------------------
mergePFASFCMAT1 <- function(dataset="heparg2d_toxcast_pfas_pe1_normal_v2",
                            file1="httr_heparg2d_toxcast_pfas_FCmat1-meanncnt0_5-plateteffect_1-shrinkage_normal.RData",
                            file2="httr_pfas_replace_heparg_FCmat1-meanncnt0_5-plateteffect_1-shrinkage_normal.RData",
                            dir="../input/fcdata/new_versions/",
                            do.load=T){
  printCurrentFunction()
  if(do.load) {
    cat("   start loading RData file\n")
    infile <- paste0(dir,file1)
    print(infile)
    load(file=infile)
    cat("finished loading\n")
    FCMAT1.1 <<- FCMAT1

    infile <- paste0(dir,file2)
    print(infile)
    load(file=infile)
    cat("finished loading\n")
    FCMAT1.2 <<- FCMAT1

  }
  cat("file 1:",nrow(FCMAT1.1),"\n")
  cat("file 2:",nrow(FCMAT1.2),"\n")
  FCMAT1 = rbind(FCMAT1.1,FCMAT1.2)
  cat("initial:",dataset,":",nrow(FCMAT1),"\n")

  temp = unique(FCMAT1[,c("chem_id","chem_name","dtxsid","pg_id","block_id")])
  names(temp) = c("sample_id","name","dtxsid","pg_id","block_id")
  file = paste0(dir,"PG_MAP_",dataset,".xlsx")
  write.xlsx(temp,file=file)
  #browser()

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
