#----------------------------------------------------------------------------------------
#' Add gene counts to the pathway catalog
#'
#' @param basedir The base directory to find the raw fold change and chemcial data
#' @param dataset The name of the datset to be used
#'
#' @return nothing is returned
#' @export
#----------------------------------------------------------------------------------------
addCountsToAssayCatalog<- function(basedir="../input/fcdata/",dataset="DMEM_6hr_pilot_normal_00") {

  file <- paste0(basedir,"FCMAT2_",dataset,".RData")
  print(file)
  load(file)
  file <- paste0(basedir,"CHEM_DICT_",dataset,".RData",sep="")
  load(file)
  rownames(CHEM_DICT) <- CHEM_DICT[,"sample_key"]

  file <- "../input/processed_pathway_data/pathway_catalog 2019-10-25.xlsx"
  pathways <- read.xlsx(file)

  file <- "../input/processed_pathway_data/PATHWAY_LIST_AllPathways.RData"
  load(file=file)

  gene.list <- colnames(FCMAT2)

  pathways$ngene <- 0
  pathways$ngene_in_data <- 0
  for(i in 1:nrow(pathways)) {
    pathway <- pathways[i,"pathway"]
    glist <- pathway_data[pathway][[1]]
    pathways[i,"ngene"] <- length(glist)
    glist <- glist[is.element(glist,gene.list)]
    pathways[i,"ngene_in_data"] <- length(glist)
  }
  file <- paste0("../input/processed_pathway_data/pathway_catalog 2019-10-31.xlsx")
  write.xlsx(pathways,file)
}
