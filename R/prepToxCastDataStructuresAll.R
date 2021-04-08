#--------------------------------------------------------------------------------------
#'
#' Create the needed data structures for just the selected chemicals and assays
#' @param date.string The data of the ToxCast data export
#' @param chems The chemical data frame
#' @param assay.list The list of assays to use
#' created the minimal data files for each of the min, max, and med points from toxboot
#--------------------------------------------------------------------------------------
prepToxCastDataStructuresAll <- function(date.string="190814",chems,assay.list) {
  printCurrentFunction()
  rownames(chems) <- chems$dtxsid
  dtxsid.list <- chems$dtxsid
  print(length(dtxsid.list))
  file <- paste0("../input/toxcast_matrix/tested_Matrix_",date.string,".csv" )
  temp <- read.csv(file,stringsAsFactors=F)
  rownames(temp) <- temp[,1]
  temp <- temp[,2:ncol(temp)]
  dtxsid.list <- dtxsid.list[is.element(dtxsid.list,rownames(temp))]
  print(length(dtxsid.list))
  temp <- temp[dtxsid.list,]
  temp[is.na(temp)] <- 0
  cs <- colSums(temp)
  mask <- cs
  mask[] <- 0
  mask[cs> (0.95*length(dtxsid.list))] <- 1
  mask[] <- 1
  temp <- temp[,mask==1]
  MAT.tested <- temp

  file <- paste0("../input/toxcast_matrix/hitc_Matrix_",date.string,".csv" )
  temp <- read.csv(file,stringsAsFactors=F)
  rownames(temp) <- temp[,1]
  temp <- temp[,2:ncol(temp)]
  temp <- temp[dtxsid.list,mask==1]
  temp[is.na(temp)] <- 0
  MAT.hitc <- temp

  file <- paste0("../input/toxcast_matrix/neglogac50_Matrix_",date.string,".csv" )
  temp <- read.csv(file,stringsAsFactors=F)
  rownames(temp) <- temp[,1]
  temp <- temp[,2:ncol(temp)]
  temp <- temp[dtxsid.list,mask==1]
  temp[is.na(temp)] <- 0
  temp[MAT.hitc==0] <- 0
  MAT.logac50 <- temp

  file <- paste0("../input/toxcast_matrix/modl_tp_Matrix_",date.string,".csv" )
  temp <- read.csv(file,stringsAsFactors=F)
  rownames(temp) <- temp[,1]
  temp <- temp[,2:ncol(temp)]
  temp <- temp[dtxsid.list,mask==1]
  temp[MAT.hitc==0] <- 0
  temp[is.na(temp)] <- 0
  MAT.top <- temp

  file <- paste0("../input/toxcast_matrix/zscore_Matrix_",date.string,".csv" )
  temp <- read.csv(file,stringsAsFactors=F)
  rownames(temp) <- temp[,1]
  temp <- temp[,2:ncol(temp)]
  temp <- temp[dtxsid.list,mask==1]
  temp[MAT.hitc==0] <- 0
  temp[is.na(temp)] <- 0
  MAT.z <- temp

  MAT.auc <- MAT.logac50 * MAT.top
  CHEMS <- chems[dtxsid.list,]

  file <- "../input/ToxCastDataPFASAll.RData"
  save(MAT.tested,MAT.hitc,MAT.logac50,MAT.top,MAT.z,MAT.auc,CHEMS,file=file)
}
