#--------------------------------------------------------------------------------------
#' Annotate the screen chemicals with refchemdb targets
#'
#'
#' @export
#--------------------------------------------------------------------------------------
annotateScreenChemicals <- function(do.load=F,
                                    basedir="../input/fcdata/",
                                    dataset="DMEM_6hr_screen_normal_pe_1"){
  printCurrentFunction()

  if(do.load) {
    file <- paste0(basedir,"CHEM_DICT_",dataset,".RData",sep="")
    print(file)
    load(file)
    rownames(CHEM_DICT) <- CHEM_DICT[,"sample_key"]
    CHEM_DICT <<- CHEM_DICT
    file <- "../input/chemicals/RefChemDB S12 Candidate reference chemicals after mode aggregation.xlsx"
    refchemdb <- read.xlsx(file)
    refchemdb <<- refchemdb
  }

  refchemdb <- refchemdb[refchemdb$support>2,]
  mat <- unique(CHEM_DICT[,c("dtxsid","casrn","name")])

  mat$target <- "-"

  refchemdb <- refchemdb[is.element(refchemdb$dsstox_substance_id,mat$dtxsid),]
  for(i in 1:nrow(mat)) {
    dtxsid <- mat[i,"dtxsid"]
    temp <- refchemdb[is.element(refchemdb$dsstox_substance_id,dtxsid),]
    if(nrow(temp)>0) {
      tlist <- sort(unique(temp$target))
      mat[i,"target"] <- paste(tlist,collapse="|")
    }
  }
  file <- "../input/chemicals/screen_chemicals_target_annoations.xlsx"
  write.xlsx(mat,file)
}
