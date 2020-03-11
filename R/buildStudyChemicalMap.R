#--------------------------------------------------------------------------------------
#' Visualize the variance by signature - this supplants signatureVariance.R
#'
#' @param min.ngene Signatures will only be saved if the number of genes is >= this value
#' @param max.ngene Signatures will only be saved if the number of genes is <= this value
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
buildStudyChemicalMap <- function(dataset="DMEM_6hr_screen_normal_pe_1"){
  printCurrentFunction()

  file <- paste0("../input/chemicals/",dataset,"_chemical_map bare.xlsx")
  mat <- read.xlsx(file)
  names(mat)[is.element(names(mat),"trt_grp_id")] <- "sample_key"
  names(mat)[is.element(names(mat),"chem_id")] <- "sample_id"
  names(mat)[is.element(names(mat),"dose_level")] <- "conc_index"

  file <- "../input/DSSTox/DSSTox_sample_map.xlsx"
  smap <- read.xlsx(file)
  smap <- smap[!is.na(smap$spid),]
  rownames(smap) <- smap$spid
  smap <- smap[mat$sample_id,]
  mat <- cbind(mat,smap[,c("dtxsid","casrn","name")])
  file <- paste0("../input/chemicals/",dataset,"_chemical_map.xlsx")
  write.xlsx(mat,file)
}
