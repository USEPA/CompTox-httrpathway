#--------------------------------------------------------------------------------------
#' buildSampleMap
#' Generate the sample_key x sample x DSSTox file
#' @param dataset Name of hte HTTr dataset
#' @param dsstox.file Name of the DSStox chemical file
#' @param dir Directory where the FCMAT1 files lives
#' @param outfile Name of the output file
#' @param do.read If TRUE, read in the input FCMAT1 file
#' @return nothing, just writes out to .xlsx
#' @importFrom openxlsx write.xlsx read.xlsx
#' @importFrom stringr str_replace_all
#' @export buildSampleMap
#'
#--------------------------------------------------------------------------------------
buildSampleMap <- function(dataset="DMEM_6hr_pilot_normal_pe_1",
                           dsstox.file="../input/DSSTox/DSSTox_sample_map.xlsx",
                           dir="../input/fcdata/",
                           outfile="../input/chemicals/HTTr_pilot_sample_map.xlsx",
                           do.read=F) {
  printCurrentFunction()
  if(do.read) {
    file <- paste0(dir,"FCMAT1_",dataset,".RDS")
    print(file)
    FCMAT1 <<- readRDS(file)
    cat("data loaded\n")
  }

  mat <- FCMAT1
  sid.list <- unique(mat$sample_key)
  temp <- str_split(sid.list,"_")
  sample.map <- as.data.frame(do.call(rbind,temp),stringsAsFactors=F)
  names(sample.map) <- c("spid","conc.index","conc.string","media","time")
  conc <- sample.map$conc.string
  units <- conc
  units <- str_replace_all(units,"\\.","")
  for(i in 0:9) units <- str_replace_all(units,as.character(i),"")

  for(unit in unique(units)) conc <- str_replace_all(conc,unit,"")
  conc <- as.numeric(conc)

  sample.map$conc <- conc
  sample.map$units <- units
  sample.map$sample_key <- sid.list
  sample.map <- sample.map[,c("sample_key","spid","conc.index","conc","units","media","time")]
  smat <- read.xlsx(dsstox.file)
  smat <- unique(smat[,c("spid","dtxsid","casrn","name")])
  smat <- smat[!is.na(smat$spid),]
  rownames(smat) <- smat$spid
  spid.list <- sample.map$spid
  missing <- spid.list[!is.element(spid.list,smat$spid)]
  if(length(missing)>0) {
    cat("missing sample IDs\n")
    print(missing)
    browser()
  }
  temp <- smat[sample.map$spid,c("dtxsid","casrn","name")]
  sample.map <- cbind(sample.map,temp)
  write.xlsx(sample.map,outfile)
}
