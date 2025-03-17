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
#-------------------------------------------------------------------------------------
extractRawCounts <- function(dir="../input/rawdata/cellatlas/") {
  printCurrentFunction()
  cwd = getwd()
  nwd = "../../httrpl/httrpl/"
  #source("../../httrpl/httrpl/Rlib/httrpl.R", chdir=T)
  #source("./Rlib/httrpl.R", chdir=T)
  #collection <- openMongo(host="fe.epa.gov", db="fe.epa.gov", username="readonly", passwd="ccte", collection="httr_mcf7_ph1")

  options(httrDefaultHost="fe.epa.gov")
  options(httrDefaultUser="readonly")
  options(httrDefaultPasswd="ccte")


  db1 = "httr_cell_atlas"
  db2 = "httr_tox21_cpp2"


  con = openMongo(db=db2,collection="httr_well_trt")

  host = "fe.epa.gov"
  db_collection = "httr_counts"
  db_collection = "httr_well_trt"
  httr_counts <- openMongo(host = host, db = db2, collection = db_collection)
  count_samples <- sort(httr_counts$distinct(key="sample_id"))


  httr_counts_qc <- openMongo(host = host, db = db2, collection = db_collection)
  count_samples <- splitIDs(DB=httr_counts_qc, ids=count_samples)
  count_samples <- unique(unlist(count_samples))
  browser()
  # httr_well is final location of the data
  #
  #DB = openMongo(db='httr_ph1',host='pb.epa.gov')
  #wells = pd.DataFrame(list(DB.htt_well.find({},dict(_id=0,probe_cnts_nrm=0,probe_cnts=0))))

  browser()
  # npg <- 48
  # for(pg in 1:npg) {
  #   cat("pg:",pg,"\n")
  #   newdir <- paste0(dir,"pg_",pg,"/")
  #   dir.create(dir, showWarnings = F)
  #   my_samples <- getWellInfo(db_host="pb.epa.gov", db_name="httr_mcf7_ph1", pg_id=pg)$sample_id
  #   browser()
  #   well_data <- getWellCounts(db_host="pb.epa.gov", db_name="httr_mcf7_ph1", sample_ids=my_samples)
  #
  #   browser()
  # }
}
