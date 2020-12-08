#source("https://bioconductor.org/biocLite.R")
#.libPaths( "C:/Program Files/R/R-3.3.3/library")
library(fastmatch)
library(reshape2)
#source("R/general_utils.R")
#options(java.parameters = "-Xmx8000m")
#library(httRws)
library(openxlsx)
#--------------------------------------------------------------------------------------
#' Prepare the raw data fro the replicated chemicals
#'
#' @param dir The directory to work in
#'
#--------------------------------------------------------------------------------------
raw.replicate.data.prep <- function(dir="input/raw_replicate_analysis/") {
  printCurrentFunction()

  file <- "input/chemicals/HTTR chemical master.xlsx"
  chems <- read.xlsx(file)
  esid.list <- sort(unique(chems[chems[,"duplicate"]==1,"epa_sample_id"]))
  rawmat <- NULL
  counter <- 0
  for(esid in esid.list) {
    counter <- counter+1
    temp <- NULL
    ptm0 <- proc.time()
    while(is.null(temp)) {
      cat("load",esid,":",counter,":",length(esid.list),"\n")
      temp <- getHTTrProbeCounts(esid)
    }
    ptm1 <- proc.time()
    print(ptm1-ptm0)
    browser()


    mat <- temp$l2fc

    factor.list <- c("chem_name","dsstox_sid","epa_sample_id","gene","probe_id","trt_grp_id")
    for(factor in factor.list) mat[,factor] <- as.character(mat[,factor])
    esid.list <- unique(mat[,"epa_sample_id"])
    mat.A <- mat[is.element(mat[,"epa_sample_id"],esid.list[1]),]
    mat.B <- mat[is.element(mat[,"epa_sample_id"],esid.list[2]),]
    sample_key <- paste(mat.A[,"epa_sample_id"],mat.A[,"conc"],mat.A[,"timeh"],sep="_")
    mat.A <- cbind(sample_key,mat.A)
    sample_key <- paste(mat.B[,"epa_sample_id"],mat.B[,"conc"],mat.B[,"timeh"],sep="_")
    mat.B <- cbind(sample_key,mat.B)

    mat.A$casrn <- unique(chems[is.element(chems[,"dsstox_sid"],dtx),"casrn"])[1]
    mat.B$casrn <- unique(chems[is.element(chems[,"dsstox_sid"],dtx),"casrn"])[1]
    fcmat1.A <- rbind(fcmat1.A,mat.A)
    fcmat1.B <- rbind(fcmat1.B,mat.B)
  }
  CHEM_DICT <- 	unique(fcmat1.A[,c("sample_key","epa_sample_id","conc","timeh","casrn","chem_name","dsstox_sid")])
  file <- paste0(dir,"CHEM_DICT.RData")
  save(CHEM_DICT,file=file)
  FCMAT1 <- fcmat1.A
  file <- paste0(dir,"FCMAT1.A.RData")
  save(FCMAT1,file=file)
  FCMAT1 <- fcmat1.B
  file <- paste0(dir,"FCMAT1.B.RData")
  save(FCMAT1,file=file)
}
