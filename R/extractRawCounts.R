#--------------------------------------------------------------------------------------
#' Get the raw counts from the Mongo database by plategroup
#'
#' THIS CODE SUPERCEDED by CODE IN HTTRPL
#'
#'
#'
#' @param dir The directory from which to read all of the raw filesatalog file
#' @param mc.cores The number of cores to use in reading the tsv files
#' @param method Either "gene" or "probe"
#' @param do.read If TRUE, read in the FCMAT1 file and place in a global.
#' @param chemical.file The required map from sample keys to chemical information
#' @param dsstox.file The information mapping chemicals to DSSTox IDs
#/
#' @return Global variables are created for the FC matrix (FCMAT2), the SE matrix (SEMAT2)
#' and the chemical dictionary (CHEM_DICT) which translates form the sample key
#' (sample_id_conc_time) to the individual components
#'
#--------------------------------------------------------------------------------------
extractRawCounts <- function(dir="../input/rawdata/mcf7_screen/") {
  printCurrentFunction()
  source("../../httrpl/httrpl/Rlib/httrpl.R", chdir=T)
  #collection <- openMongo(host="fe.epa.gov", db="fe.epa.gov", username="readonly", passwd="ccte", collection="httr_mcf7_ph1")

  options(httrDefaultHost="fe.epa.gov")
  options(httrDefaultUser="readonly")
  options(httrDefaultPasswd="ccte")


  npg <- 48
  for(pg in 1:npg) {
    cat("pg:",pg,"\n")
    newdir <- paste0(dir,"pg_",pg,"/")
    dir.create(dir, showWarnings = F)
    my_samples <- getWellInfo(db_host="fe.epa.gov", db_name="httr_mcf7_ph1", pg_id=pg)$sample_id
    well_data <- getWellCounts(db_host="fe.epa.gov", db_name="httr_mcf7_ph1", sample_ids=my_samples)

    browser()



  }

}
