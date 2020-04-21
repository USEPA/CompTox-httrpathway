#--------------------------------------------------------------------------------------
#' Build the FCMAT1 data set
#' This is a variant on the original buildFCMAT1 to handle large data sets
#' hopefully a faster version of the very for-loopy FCMAT1 found in the myGSEA pagackge taken from Derik Hggard
#' on feb 20 2020.  Written by Bryant Chambers mar 15 2020
#'
#' about:
#' this script builds the giant FCMAT1 from a folder of tsv files in a directory of your choosing
#'
#' important. it asssumes that you have all files as .tsv s in a single folder. it doesn't check
#' that all files are tsvs or exclude non tsv files. I just wrote this as a trial and it seems to
#' significantly outperform the original. (a 3421 file set took 6 min while the original took 8 hours)
#'
#' buildFCMAT1.v3(dataset="DMEM_6hr_screen_normal_pe_1",
#' dir="../input/httr_mcf7_screen/meanncnt0_5-plateteffect_0-shrinkage_normal_DMEM_6/")
#'
#' @param dataset The name to give to the data set
#' @param dir The directory from which to read all of the raw filesatalog file
#' @param mc.cores The number of cores to use in reading the tsv files
#' @return A file with the FCMAT1 data is written to "../input/fcdata/FCMAT1_",dataset,".RData"
#'
#' @export
library(tidyverse)
#--------------------------------------------------------------------------------------
buildFCMAT1.v3 <- function(dataset="DMEM_6hr_screen_normal_pe_1",
                          dir="../input/httr_mcf7_screen/meanncnt0_5-plateteffect_1-shrinkage_normal_DMEM_6/",
                          mc.cores=1){
#  buildFCMAT1.v3 = function(dataset="DMEM_6hr_pilot_normal_pe_1",
#                            dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_1-shrinkage_normal/DMEM_6/",
#                            mc.cores=1){
    printCurrentFunction()

  filelist <- list.files(dir) #list of all files to iterate through
  ###########################################################################################################
  #read all files and bind them simultaneously
  #this builds a perliminary version of the FCMAT1
#  allframes <- suppressMessages(lapply(1:length(filelist),function(x)read_tsv(paste0(dir,filelist[x]))))
  allframes <- suppressMessages(mclapply(1:length(filelist),function(x)read_tsv(paste0(dir,filelist[x])),mc.cores=mc.cores))
  cat("   finish reading in all files:",length(allframes),":",length(filelist),"\n")
  if(length(allframes)!=length(allframes)) {
    cat("*** Not all data read in\n")
    browser()
  }
  #allframesbound <- do.call(rbind, allframes)
  allframesbound <- rbindlist(allframes)
  cat("   finish rbind:",dim(allframesbound),"\n")

  #allframesbound <- unique(allframesbound)
  #cat("   unique data:",dim(allframesbound),"\n")

  ###########################################################################################################
  #clean up the mat so it look identical to the original fcmat
  names(allframesbound) <- c("probe_id","sample_key","basemean","l2fc","se","stat","pvalue","padj")               ## double check these so that they match yours
  allframesbound$basemean <- as.integer(allframesbound$basemean)
  FCMAT1 <- as.data.frame(allframesbound)
  cat("   finish builing final data frame\n")
  file <- paste0("../input/fcdata/FCMAT1_",dataset,".v3.RData")
  cat("   start saving RData file\n")
  save(FCMAT1,file=file)
  cat("finish saving RData file\n")
}
