#--------------------------------------------------------------------------------------
#' get the mapping between the plate groups and the samples
#'
#' version to start with Logan's database export
#' The difference between this version and the original is that there are extra columns
#' The function just changes one column name and writes the file to a standard name and place
#'
#' @param dataset The name to give to the data set
#' @param dir The directory from which to read all of the raw files
#' @param do.load if T, load the initial file
#' @return A mapping the sampels to the plate groups
#' f
#'
#' @export
#--------------------------------------------------------------------------------------
pg_id.to.sample_id <- function(do.load=F,
                               dataset="mcf7_ph1_pe1_normal_block_123_allPG",
                               dir="../input/fcdata/"){
  printCurrentFunction()
  if(do.load) {
    cat("   start loading RData file\n")
    infile <- paste0(dir,"FCMAT1_",dataset,".RData")
    print(infile)
    load(file=infile)
    cat("finished loading\n")
    FCMAT1 <<- FCMAT1
  }
  fcmat1 = FCMAT1

  temp = unique(fcmat1[,c("chem_id","chem_name","dtxsid","pg_id","block_id")])
  names(temp) = c("sample_id","name","dtxsid","pg_id","block_id")

  file = paste0(dir,"PG_MAP_",dataset,".xlsx")
  write.xlsx(temp,file=file)
}
