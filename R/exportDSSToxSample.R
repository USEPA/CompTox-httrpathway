#--------------------------------------------------------------------------------------
#'
#' Generate the sample x DSSTox file
#'
#' @param outfile Name of the file to be written
#'
#--------------------------------------------------------------------------------------
exportDSSToxSample <- function(outfile="../input/DSSTox/DSSTox_sample_map.xlsx") {
  printCurrentFunction()

  if(!exists("DSSTOX")) {
    file <- "../input/DSSTox/DSSTox 2020-01-02.RData"
    load(file=file)
    DSSTOX <<- DSSTOX
  }
  dsstox <- DSSTOX
  rownames(dsstox) <- dsstox$dtxsid

  query <- paste0("
  SELECT
    scd.blinded_sample_id AS spid, gs.dsstox_substance_id as dtxsid
  FROM
    ro_prod_dsstox.generic_substances gs
        JOIN
    ro_prod_chemtrack.bottles b ON b.efk_generic_substance_id = gs.id
        JOIN
    ro_prod_chemtrack.shipped_chemical_details scd ON scd.fk_bottle_id = b.id
  ")
  res <- runQuery(query,"ro_prod_dsstox")
  res <- unique(res)
  dsstox <- dsstox[res$dtxsid,c("casrn","name")]
  res <- cbind(res,dsstox)
  write.xlsx(res,outfile)
}
