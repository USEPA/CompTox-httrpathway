#--------------------------------------------------------------------------------------
#'
#' Build lane plots by chemical list and signature class, across the datasets
#' @param to.file If TRUE, write plots to a file
#--------------------------------------------------------------------------------------
signaturePOD <- function(sigset="pilot_large_all_100CMAP",
                       dataset="DMEM_6hr_pilot_normal_pe_1",
                       method="mygsea",
                       hccut=0.5,
                       do.load=F) {
  printCurrentFunction()

  if(do.load) {
    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat <- SIGNATURE_CR
    MAT <<- mat
  }
  mat = MAT
  result <- unique(mat[,c("dtxsid","casrn","name")])
  result <- result[order(result$name),]
  result$nsignature <- NA
  result$signature_pod_min <- NA
  result$signature_pod_min.lci <- NA
  result$signature_pod_min.uci <- NA
  result$signature_pod_95 <- NA
  result$signature_pod_95.lci <- NA
  result$signature_pod_95.uci <- NA
  result$signature_pod_95.count <- NA
  result$signature_pod_95.count.lci <- NA
  result$signature_pod_95.count.uci <- NA
  rownames(result) <- result$dtxsid
  for(dtxsid in result$dtxsid) {
    temp <- mat[is.element(mat$dtxsid,dtxsid),]
    temp <- temp[temp$hitcall>hccut,]
    npath <- nrow(temp)
    result[dtxsid,"nsignature"] <- npath
    temp <- temp[order(temp$bmd),]
    bmdl <- temp$bmdl
    bmdu <- temp$bmdu
    ratio <- bmdu/bmdl
    ratio [is.na(ratio)] <- 100
    temp <- temp[ratio<40,]
    result[dtxsid,"signature_pod_min"] <- temp[1,"bmd"]
    result[dtxsid,"signature_pod_min.lci"] <- temp[1,"bmdl"]
    result[dtxsid,"signature_pod_min.uci"] <- temp[1,"bmdu"]
    minval <- 5
    result[dtxsid,"signature_pod_95.count"] <- temp[minval,"bmd"]
    result[dtxsid,"signature_pod_95.count.lci"] <- temp[minval,"bmdl"]
    result[dtxsid,"signature_pod_95.count.uci"] <- temp[minval,"bmdu"]
    if(npath>100) minval <- 0.05*npath
    result[dtxsid,"signature_pod_95"] <- temp[minval,"bmd"]
    result[dtxsid,"signature_pod_95.lci"] <- temp[minval,"bmdl"]
    result[dtxsid,"signature_pod_95.uci"] <- temp[minval,"bmdu"]
  }
  x <- result$signature_pod_min
  x[is.na(x)] <- 1000
  result$signature_pod_min <- x

  x <- result$signature_pod_95
  x[is.na(x)] <- 1000
  result$signature_pod_95 <- x

  x <- result$signature_pod_95.count
  x[is.na(x)] <- 1000
  result$signature_pod_95.count <- x

  file <- paste0("../output/signature_pod/signature_pod_",sigset,"_",dataset,"_",method,"_",hccut,".xlsx")
  write.xlsx(result,file)
}

