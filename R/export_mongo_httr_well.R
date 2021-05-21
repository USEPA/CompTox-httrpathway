library(mongolite)
source("./httrpl/Rlib/httrpl.R",chdir=T)
#--------------------------------------------------------------------------------------
#' Get the raw counts from the Mongo database
#'
#' @param db The name of the Mongo database
#' @param collection THe name of the collection to export
#' @param dir The directory where the data will be stored
#'
#' Collections
#' * httr_cell_atlas
#' * httr_tox21_cpp2
#' @export
#-------------------------------------------------------------------------------------
export_mongo_httr_well <- function(db="httr_cell_atlas",collection="httr_well_trt",dir="../input/rawdata/cellatlas/") {
  printCurrentFunction()

  options(httrDefaultHost="fe.epa.gov")
  options(httrDefaultUser="readonly")
  options(httrDefaultPasswd="ccte")
  #username = "readonly"
  #password = "ccte"
  my_host="fe.epa.gov"

  #collection <- openMongo(host=host,db=db,user=username,passwd=password,collection=db)
  collection <- openMongo(db=db,collection=collection)
  mat1 = getWellInfo(db_name=db,db_host=my_host)
  sid.list = mat1$sample_id
  welldata = getWellCounts(db_host=my_host, db_name=db, sample_ids=sid.list)
  info = welldata[[1]]
  counts = welldata[[2]]
  for(i in 1:ncol(info)) if(class(info[,i])=="factor") info[,i] = as.character(info[,i])

  file = paste0(dir,db,"_info.RData")
  save(info,file=file)
  file = paste0(dir,db,"_counts.RData")
  save(counts,file=file)
}
