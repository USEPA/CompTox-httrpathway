#--------------------------------------------------------------------------------------
#' Build the chemical mapping file for the MCF7 screen
#'
#' @param dataset The name to give to the data set
#' @param dir The directory from which to read all of the raw filesatalog file
#' @param filetype Either tsv or RData
#' @return A file with the FCMAT1 data is written to "../input/fcdata/FCMAT0_",dataset,".RData"
#' @export
#--------------------------------------------------------------------------------------
buildChemicalMappingFile = function(dir="../input/chemicals/"){
  printCurrentFunction()
  file <- paste0(dir,"httr_mcf7_ph1_chem.xlsx")
  mat1 <- read.xlsx(file)
  file <- paste0(dir,"HTTr.MCF7.Screen.Sample.Key.20180605.xlsx")
  mat2 <- read.xlsx(file)

  mat2$chem_name <- NA
  mat2$casrn <- NA
  mat2$dtxsid <- NA





  browser()
}
