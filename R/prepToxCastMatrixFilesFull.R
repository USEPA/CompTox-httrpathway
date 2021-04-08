library(tcpl)
library(openxlsx)
library(DBI)
library(RMySQL)
library(stringi)
library(stringr)
#--------------------------------------------------------------------------------------
#'
#' Prepare the full ToxCast matrix files
#' Input files are ...
#' chemicals: "../input/PFAS 149 2019-08-14.xlsx"
#' assays: "../input/assay_list.txt"
#' toxboot information: "../input/pfas_toxboot.RData"
#' The raw toxcast matrices are put into "../input/toxcast_matrix_full/"
#'
#' All concentrations are in log10(uM)
#' @param invitrodb.extract If TRUE, extrac hte raw data from invitrodb, else
#' created the minimal data files for each of the min, max, and med points from toxboot
#' @param invitrodb.extract.flags If TRUE, extract the flag data
#' @param invitrodb.extract.cytotox If TRUE, extract the cytoxicity data
#' @param prep.toxboot If TRUE, extract the toxboot data
#--------------------------------------------------------------------------------------
prepToxCastMatrixFilesFull <- function(invitrodb.extract=F,
                                       invitrodb.extract.flags=F,
                                       invitrodb.extract.cytotox=F,
                                       prep.toxboot=F) {
  printCurrentFunction()
  tcplConf(user='rjudson', pass='Catman2@', host='ccte-mysql-res.epa.gov', drvr = 'MySQL',db = 'invitrodb')

  odir <- "../input/toxcast_matrix/"

  if(invitrodb.extract) {
    #tmp1 <- tcplLoadAeid(fld="internal_ready",val=1)
    tmp1 <- tcplLoadAeid()
    tmp1 <- as.data.frame(tmp1)
    #tmp1 <- tmp1[is.element(tmp1[,"aenm"],assay.list),]
    aeid.list <- tmp1[,"aeid"]

    tmp1 <- tcplLoadChem()
    tmp1 <- as.data.frame(tmp1)
    tmp1 <- unique(tmp1[,c("dsstox_substance_id","chid")])
    #tmp1 <- tmp1[is.element(tmp1$dsstox_substance_id,chems$dtxsid),]
    chid.list <- tmp1$chid

    temp <- tcplVarMat(chid = chid.list, aeid = aeid.list, row.id = "dsstox_substance_id",
                       flag = TRUE, cyto.pars = list(), include.na.chid = FALSE, odir = odir,file.prefix = NULL,
                       add.vars=c("resp_max","resp_min","max_mean","max_mean_conc","max_med","max_med_conc","logc_max","logc_min",
                                  "hill_tp","hill_ga","hill_gw","hill_er","gnls_tp","gnls_ga","gnls_gw","gnls_la","gnls_lw",
                                  "modl_tp","modl_ga","modl_gw","modl_la","modl_lw","modl_acc","modl_acb","modl_ac10"))
  }
  if(invitrodb.extract.flags) {
    mat <- tcplPrepOtpt(tcplLoadData(6L))
    file <- paste(odir,"toxcast_flags_",Sys.Date(),".csv",sep="")
    write.csv(mat,file)
  }
  if(invitrodb.extract.cytotox) {
    mat <- runQuery(
      "select a.dsstox_substance_id as dtxsid, a.casn as casrn, a.chnm as name,
      b.cytotox_median_raw,b.cytotox_mad,b.global_mad,b.cytotox_median_log,b.cytotox_median_um,b.cytotox_lower_bound_um,b.ntested,b.nhit,b.cytotox_lower_bound_log
      from chemical a, cytotox b
      where a.chid=b.chid","invitrodb"
    )
    file <- paste(odir,"cytotox_ranges_",Sys.Date(),".xlsx",sep="")
    write.xlsx(mat,file)
  }
  if(prep.toxboot) {
    prepToxBootFile(chems,assay.list)
  }
}
