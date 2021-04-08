library(tcpl)
library(openxlsx)
library(DBI)
library(RMySQL)
library(stringi)
library(stringr)
#--------------------------------------------------------------------------------------
#'
#' Export the concentration response data for selected assays
#' The results are put into "../input/toxcast_matrix/"
#'
#' All concentrations are in log10(uM)
#' @param assay.list The list of assays to use
#' @param invitrodb.extract If TRUE, extrac hte raw data from invitrodb, else
#' created the minimal data files for each of the min, max, and med points from toxboot
#' @param invitrodb.extract.flags If TRUE, extract the flag data
#' @param invitrodb.extract.cytotox If TRUE, extract the cytoxicity data
#' @param prep.toxboot If TRUE, extract the toxboot data
#--------------------------------------------------------------------------------------
prepToxCastConcResp <- function(infile="../input/toxcast_matrix/ppar_assays.xlsx",
                                outfile="../input/toxcast_matrix/ppar_cr.xlsx") {
  printCurrentFunction()
  tcplConf(user='rjudson', pass='Catman2@', host='ccte-mysql-res.epa.gov', drvr = 'MySQL',db = 'invitrodb')

  dir <- "../input/toxcast_matrix/"
  assays = read.xlsx(infile)
  tmp1 <- tcplLoadAeid()
  tmp1 <- as.data.frame(tmp1)
  tmp1 <- tmp1[is.element(tmp1[,"aenm"],assays$aenm),]
  aeid.list <- tmp1[,"aeid"]
  #aeid.list = aeid.list[is.element(aeid.list,assays$aeid)]

  #tmp1 <- tcplLoadChem()
  #tmp1 <- as.data.frame(tmp1)
  #tmp1 <- unique(tmp1[,c("dsstox_substance_id","chid")])
  #tmp1 <- tmp1[is.element(tmp1$dsstox_substance_id,chems$dtxsid),]
  #chid.list <- tmp1$chid
  mat = tcplLoadData(lvl=4,fld="aeid",val=aeid.list,type="mc")
  mat = as.data.frame(mat)
  browser()
}
