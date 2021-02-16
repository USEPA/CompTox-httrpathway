#--------------------------------------------------------------------------------------
#' pull our HTTr data for teh Unilever chemicals
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#'
#--------------------------------------------------------------------------------------
Unilever.httr.merge <- function() {
  printCurrentFunction()
  celltype = "HepaRG"
  file = paste0("../input/Unilever/Unilever HTTr overlap ",celltype,".xlsx")
  mat = read.xlsx(file)
  for(i in 5:7) names(mat)[i] = paste0(celltype,".",names(mat)[i])
  heparg = mat

  celltype = "MCF7"
  file = paste0("../input/Unilever/Unilever HTTr overlap ",celltype,".xlsx")
  mat = read.xlsx(file)
  for(i in 5:7) names(mat)[i] = paste0(celltype,".",names(mat)[i])
  mcf7 = mat

  celltype = "U2OS"
  file = paste0("../input/Unilever/Unilever HTTr overlap ",celltype,".xlsx")
  mat = read.xlsx(file)
  for(i in 5:7) names(mat)[i] = paste0(celltype,".",names(mat)[i])
  u2os = mat

  res = cbind(heparg[,c(1,2,3,5:7)],mcf7[,c(5:7)],u2os[,c(5:7)])

  file = paste0("../input/Unilever/Unilever HTTr overlap merged.xlsx")
  write.xlsx(res,file)


  file = "../input/chemicals/HTTR_chemical_annotations 2020-12-02.xlsx"
  chems = read.xlsx(file)
  chems = chems[is.element(chems$dtxsid,res$dtxsid),]
  file = paste0("../input/Unilever/Unilever use information.xlsx")
  write.xlsx(chems,file)
}

