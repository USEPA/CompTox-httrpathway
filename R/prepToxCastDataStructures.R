#--------------------------------------------------------------------------------------
#'
#' Create the needed data structures for just the selected chemicals and assays
#' @param date.string The data of the ToxCast data export
#' @param chems The chemical data frame
#' @param assay.list The list of assays to use
#' created the minimal data files for each of the min, max, and med points from toxboot
#--------------------------------------------------------------------------------------
prepToxCastDataStructures <- function(date.string="210310") {
  printCurrentFunction()

  file = "../input/toxcast_matrix/assay_info_2020-12-17.xlsx"
  assay.data = read.xlsx(file)
  source.list = c("ACEA","ACEA_Cytotoxicity","APR_Cytotoxicity","APR_dn","APR_up",
                  "ATG_CIS","ATG_TRANS","BSK_Cytotoxicity","BSK_down","BSK_up","ZF_terata")
  assay.list = assay.data[is.element(assay.data$source,source.list),"assay"]

  file <- paste0("../input/toxcast_matrix/tested_Matrix_",date.string,".csv" )
  temp <- read.csv(file,stringsAsFactors=F)
  rownames(temp) <- temp[,1]
  temp <- temp[,2:ncol(temp)]
  dtxsid.list = rownames(temp)
  temp <- temp[dtxsid.list,]
  assay.list = assay.list[is.element(assay.list,names(temp))]
  temp = temp[,assay.list]
  MAT.tested <- temp

  file <- paste0("../input/toxcast_matrix/hitc_Matrix_",date.string,".csv" )
  temp <- read.csv(file,stringsAsFactors=F)
  rownames(temp) <- temp[,1]
  temp <- temp[,2:ncol(temp)]
  temp[is.na(temp)] <- 0
  temp = temp[,assay.list]
  MAT.hitc <- temp

  file <- paste0("../input/toxcast_matrix/neglogac50_Matrix_",date.string,".csv" )
  temp <- read.csv(file,stringsAsFactors=F)
  rownames(temp) <- temp[,1]
  temp <- temp[,2:ncol(temp)]
  temp[is.na(temp)] <- 0

  temp = temp[,assay.list]
  temp[MAT.hitc==0] <- 0
  temp[MAT.tested==0] <- 0
  MAT.logac50 <- temp
  x = MAT.logac50
  x=cbind(rownames(x),x)
  names(x)[1] = "dtxsid"
  y = reshape2::melt(x,id.vars="dtxsid",value.name="logac50")
  y = y[y$logac50>0,]
  val = 6-y$logac50
  val = 10**val
  y$logac50 = val
  names(y)[2] = "assay"
  names(y)[3] = "ac50"
  y$source = "-"
  for(source in source.list) {
    assay.list = assay.data[is.element(assay.data$source_group,source),"assay"]
    y[is.element(y$assay,assay.list),"source"] = source
  }

  TOXCAST = y
  file = "../input/toxcast_matrix/toxcast_active_by_source.RData"
  save(TOXCAST,file=file)
}
