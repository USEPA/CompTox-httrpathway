#--------------------------------------------------------------------------------------
#' Extract annotations for CMAP chemicals
#'
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
CMAP_RefchemDB <- function(){
  printCurrentFunction()

  file <- "../input/signatures/CMAP refchemdb input.xlsx"
  input <- read.xlsx(file)
  input$dtxsid <- NA
  for(i in 1:nrow(input)) {
    x <- input[i,1]
    if(contains(x,"DTX")) {
      y <- str_split(x," ")[[1]]
      for(j in 1:length(y)) {
        if(substr(y[j],1,3)=="DTX") {
          input[i,"dtxsid"] <- y[j]
        }
      }
    }
  }

  dtxsid.list <- unique(input$dtxsid)
  file <- "../../RefChemDB/gene_target_info/chemical_gene_target_info_combined_2017-03-03_manual.xlsx"
  refchem <- read.xlsx(file)

  if(!exists("DSSTOX"))  {
    file <- "../input/DSSTox/DSSTox 2020-01-02.RData"
    load(file=file)
    DSSTOX <<- DSSTOX
  }
  dsstox <- DSSTOX
  dsstox <- dsstox[is.element(dsstox$dtxsid,dtxsid.list),]
  rownames(dsstox) <- dsstox$casrn

  casrn.list <- refchem$casrn
  dt <- dsstox[casrn.list,c("dtxsid","casrn","name")]
  refchem$dtxsid <- dt$dtxsid

  input$target <- NA
  for(i in 1:nrow(input)) {
    dtxsid <- input[i,"dtxsid"]
    if(!is.na(dtxsid)) {
      if(is.element(dtxsid,refchem$dtxsid)) {
        temp <- refchem[is.element(refchem$dtxsid,dtxsid),]
        target <- paste(temp$gene_symbol_original, collapse="|")
        input[i,"target"] <- target
      }
    }
  }
  file <- "../input/signatures/CMAP refchemdb output.xlsx"
  write.xlsx(input,file)
  browser()
}
