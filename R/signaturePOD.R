#--------------------------------------------------------------------------------------
#'
#' Calculate PODs at the signature level
#' @param do.laod If TRUE, load the input data into memory
#' @param sigset Name of signature set.
#' @param dataset Name of data set.
#' @param method Pathway scoring method in c("fc", "gsva", "gsea")
#' @param bmr_scale	bmr scaling factor. Default = 1.349
#' @param hccut Remove rows with hitcall less than this value
#'
#' @export
#--------------------------------------------------------------------------------------
signaturePOD <- function(do.load=F,
                         sigset="screen_large",
                         dataset="PFAS_U2OS",
                         method="fc",
                         bmr_scale=1.349,
                         hccut=0.95) {
  printCurrentFunction()

  if(do.load) {
    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    if(bmr_scale!=1.349) file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_bmr_scale_",bmr_scale,"_0.05_conthits.RData")
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
    if(nrow(temp)>1) {
      minval <- 5
      if(minval>nrow(temp)) minval = nrow(temp)
      result[dtxsid,"signature_pod_95.count"] <- temp[minval,"bmd"]
      result[dtxsid,"signature_pod_95.count.lci"] <- temp[minval,"bmdl"]
      result[dtxsid,"signature_pod_95.count.uci"] <- temp[minval,"bmdu"]
      minval <- 0.05*npath
      if(minval<5) minval=5
      if(minval>nrow(temp)) minval = nrow(temp)
      result[dtxsid,"signature_pod_95"] <- temp[minval,"bmd"]
      result[dtxsid,"signature_pod_95.lci"] <- temp[minval,"bmdl"]
      result[dtxsid,"signature_pod_95.uci"] <- temp[minval,"bmdu"]
    }
  }

  #result <- result[!is.na(result$signature_pod_95.count),]
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
  if(bmr_scale!=1.349) file <- paste0("../output/signature_pod/signature_pod_",sigset,"_",dataset,"_",method,"_bmr_scale_",bmr_scale,"_",hccut,".xlsx")
  write.xlsx(result,file)
}

