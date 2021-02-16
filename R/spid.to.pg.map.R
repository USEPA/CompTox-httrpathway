#--------------------------------------------------------------------------------------
#' save the mapping from plate group to sample_id
#'
#'
#' @param dataset The name to give to the data set
#' @param dir The directory from which to read all of the raw files
#' @param infile The nae of the input file
#' @return A file with the teh unique sampel_id to pg_id map written to "../input/fcdata/spid.to.pg.map ",dataset,".xlsx"
#'
#' @export
library(tidyverse)
#--------------------------------------------------------------------------------------
spid.to.pg.map <- function(dataset="mcf7_ph1_pe1_normal_block_123",
                           dir="../input/fcdata/new_versions/",
                           infile="httr_mcf7_ph1_FCmat1_meanncnt0_5-plateteffect_1-shrinkage_normal_good_pg.RData",
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
  mat = unique(FCMAT1[,c("chem_id","pg_id","dtxsid","chem_name")])
  names(mat) = c("spid","pg_id","dtxsid","name")
  file = paste0("../input/fcdata/spid.to.pg.map ",dataset,".xlsx")
  write.xlsx(mat,file)
  browser()
}
