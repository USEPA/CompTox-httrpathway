library(openxlsx)
library(tidyr)
library(stringr)
library(stringi)
library(reshape2)
library(gplots)
library(RColorBrewer)
#--------------------------------------------------------------------------------------
#' Prepare the input file for the refchemDB analyses
#--------------------------------------------------------------------------------------
refchemdbPrep <- function(do.load=F) {
  printCurrentFunction()
  if(do.load) {
    file = "../input/chemicals/RefChemDB S7 The filtered contents of the database.xlsx"
    file <- "../input/chemicals/S11 Candidate reference chemicals after mode aggregation.xlsx"
    RCDB <<- read.xlsx(file)
  }
  dataset.list = c(
    "heparg2d_toxcast_pfas_pe1_normal",
    "mcf7_ph1_pe1_normal_good_pg",
    "u2os_toxcast_pfas_pe1_normal"
  )
  chems <- NULL
  for(dataset in dataset.list) {
    file = paste0("../input/fcdata/CHEM_DICT_",dataset,".RData")
    load(file=file)
    chems = rbind(chems,CHEM_DICT)
  }
  chems = unique(chems[,c("dtxsid","name","casrn")])
  rcdb = RCDB
  rcdb <- rcdb[rcdb$count>=3,]
  rcdb = rcdb[is.element(rcdb$dsstox_substance_id,chems$dtxsid),]
  names(rcdb)[1] = "dtxsid"
  file = "../input/signatures/refchemdb_chem_filtered.xlsx"
  write.xlsx(rcdb,file)
}
