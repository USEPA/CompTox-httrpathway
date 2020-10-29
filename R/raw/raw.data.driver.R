library(fastmatch)
library(reshape2)
#library(httRws)
library(openxlsx)
library(tidyr)
library(grDevices)
library(RColorBrewer)
library(gplots)
library(ape)
#source("R/general_utils.R")
#options(java.parameters = "-Xmx8000m")
#--------------------------------------------------------------------------------------
#' run all of the raw data QC aanlyses
#'
#' @param dir The directory to work in
#'
#--------------------------------------------------------------------------------------
raw.data.driver <- function(do.data.prep=F,
                            do.data.analysis.1=F,
                            do.data.analysis.2=F,
                            do.data.analysis.3=F,
                            do.data.analysis.4=F,
                            dirin="../../../Raw/ctrl/",
                            dirout="analysis/raw_count_analysis/data/",
                            cutoff=0.95) {
  printCurrentFunction()
  if(do.data.prep) {
    welltype.list <- c("bl_dmso","bl_tsa","hbrr","uhrr","untreated","vehicle_control","refchems")
    for(welltype in welltype.list) raw.data.prep(welltype=welltype,dirin=dirin,dirout=dirout,cutoff=cutoff)
  }

  if(do.data.analysis.1) {
    raw.data.analysis(to.file=T,welltype="vehicle_control",do.read.dmso=T,do.read.other=F,prep.l2fcmat=F,
                      analysis_id=1,cutoff=cutoff,l2fcmax=10,
                      dirin="analysis/raw_count_analysis/data/",
                      dirout="analysis/raw_count_analysis/")
  }
  if(do.data.analysis.2) {
    welltype <- "bl_dmso"
    raw.data.analysis(to.file=T,welltype=welltype,do.read.dmso=T,do.read.other=T,prep.l2fcmat=T,
                      analysis_id=2,
                      cutoff=cutoff,l2fcmax=10,read_id=-1,
                      dirin="analysis/raw_count_analysis/data/",
                      dirout="analysis/raw_count_analysis/")
    welltype.list <- c("bl_tsa","hbrr","uhrr","untreated","GEN","TSA","SIRO")
    for(welltype in welltype.list) {
      raw.data.analysis(to.file=T,welltype=welltype,do.read.dmso=F,do.read.other=T,prep.l2fcmat=T,
                        analysis_id=2,
                        cutoff=cutoff,l2fcmax=10,read_id=-1,
                        dirin="analysis/raw_count_analysis/data/",
                        dirout="analysis/raw_count_analysis/")
    }
  }
  if(do.data.analysis.3) {
    raw.data.analysis(to.file=T,welltype="bl_tsa",do.read.dmso=F,do.read.other=F,prep.l2fcmat=F,
                      analysis_id=3,
                      read_id=3,welltypeA="hbrr",welltypeB="uhrr",
                      cutoff=cutoff,l2fcmax=10,
                      dirin="analysis/raw_count_analysis/data/",
                      dirout="analysis/raw_count_analysis/")
    raw.data.analysis(to.file=T,welltype="bl_tsa",do.read.dmso=F,do.read.other=F,prep.l2fcmat=F,
                      analysis_id=3,
                      read_id=3,welltypeA="bl_dmso",welltypeB="bl_tsa",
                      cutoff=cutoff,l2fcmax=10,
                      dirin="analysis/raw_count_analysis/data/",
                      dirout="analysis/raw_count_analysis/")

  }
  if(do.data.analysis.4) {
    raw.data.analysis(to.file=T,welltype="bl_tsa",do.read.dmso=F,do.read.other=F,prep.l2fcmat=F,
                      analysis_id=4,
                      cutoff=cutoff,l2fcmax=10,
                      dirin="analysis/raw_count_analysis/data/",
                      dirout="analysis/raw_count_analysis/")
  }
}
