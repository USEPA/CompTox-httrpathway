#--------------------------------------------------------------------------------------
#'
#' Build lane plots by chemical list and pathway class, across the datasets
#' @param to.file If TRUE, write plots to a file
#--------------------------------------------------------------------------------------
pathwayPOD <- function(pathset="PathwaySet_20191107",
                       dataset="DMEM_6hr_pilot_normal_pe_0",
                       method = "fc",
                       hit.threshold=0.5) {
  printCurrentFunction()

  file <- paste0("../output/pathway_conc_resp_summary/PATHWAY_CR_",pathset,"_",dataset,"_",method,"_0.05_conthits.RData")
  print(file)
  load(file=file)
  mat <- PATHWAY_CR

  result <- unique(mat[,c("dtxsid","casrn","name")])
  result <- result[order(result$name),]
  result$pathway_pod_min <- NA
  result$pathway_pod_95 <- NA
  rownames(result) <- result$dtxsid
  for(dtxsid in result$dtxsid) {
    temp <- mat[is.element(mat$dtxsid,dtxsid),]
    temp <- temp[temp$hitcall>hit.threshold,]
    npath <- nrow(temp)
    temp <- temp[order(temp$bmd10),]
    result[dtxsid,"pathway_pod_min"] <- temp[1,"bmd10"]
    result[dtxsid,"pathway_pod_95"] <- temp[0.05*npath+1,"bmd10"]
  }
  x <- result$pathway_pod_min
  x[is.na(x)] <- 1000
  result$pathway_pod_min <- x

  x <- result$pathway_pod_95
  x[is.na(x)] <- 1000
  result$pathway_pod_95 <- x

  file <- paste0("../output/pathway_pod/pathway_pod_",pathset,"_",dataset,"_",method,".xlsx")
  write.xlsx(result,file)
}

