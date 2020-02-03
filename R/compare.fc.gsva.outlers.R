#--------------------------------------------------------------------------------------
#'
#' Compare the outliers from the FC and GSVA calculations
#' @param to.file If TRUE, write plots to a file
#--------------------------------------------------------------------------------------
compare.fc.gsva.outlers <- function(to.file=F,
                                    method2="gsva",
                                    dataset="DMEM_6hr_pilot_none_pe_1",
                                    pathset="PathwaySet_20191107") {
  printCurrentFunction()
  if(to.file) {
    fname <- paste0("../output/pod_fc_gsva_compare/pod_fc_gsva_compare_outliers_",dataset,"_",pathset,"_",method2,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,2,2))
  file <- "../input/cytotoxicity summary wide allchems.xlsx"
  CYTOTOX <- read.xlsx(file)
  rownames(CYTOTOX) <- CYTOTOX$dtxsid

  method <- "fc"
  file <- paste0("../output/pathway_conc_resp_summary/PATHWAY_CR_",pathset,"_",dataset,"_",method,"_0.05_conthits.RData")
  print(file)
  load(file)
  fc <- PATHWAY_CR
  method <- method2
  file <- paste0("../output/pathway_conc_resp_summary/PATHWAY_CR_",pathset,"_",dataset,"_",method,"_0.05_conthits.RData")
  print(file)
  load(file)
  gsva <- PATHWAY_CR

  rownames(fc) <- paste(fc$dtxsid,fc$pathway)
  rownames(gsva) <- paste(gsva$dtxsid,fc$pathway)

  temp1 <- fc[fc$hitcall>0.5,]
  temp1 <- temp1[temp1$bmd<1,]
  temp2 <- gsva[gsva$hitcall>0.5,]
  temp2 <- temp2[temp2$bmd<1,]

  rn1 <- rownames(temp1)
  rn2 <- rownames(temp2)
  rn <- unique(c(rn1,rn2))
  fc <- fc[rn,]
  gsva <- gsva[rn,]

  mat <- fc[,c("dtxsid","casrn","name","pathway")]
  mat$bmd.fc <- fc$bmd
  mat$hitcall.fc <- fc$hitcall
  mat$bmd.gsva <- gsva$bmd
  mat$hitcall.gsva <- gsva$hitcall
  x <- mat$bmd.fc
  x[is.na(x)] <- 1000
  mat$bmd.fc <- x
  x <- mat$bmd.gsva
  x[is.na(x)] <- 1000
  mat$bmd.gsva <- x

  mat$delta <- abs(mat$bmd.fc-mat$bmd.gsva)

  file <- paste0("../output/pod_fc_gsva_compare/pod_fc_gsva_compare_outliers_",dataset,"_",pathset,"_",method2,".xlsx")
  write.xlsx(mat,file)

  do.plot.1 <- T
  if(do.plot.1) {
    dtxsid.list <- unique(fc$dtxsid)
    for(dtxsid in dtxsid.list) {
      temp <- mat[is.element(mat$dtxsid,dtxsid),]
      rn <- rownames(temp)
      fc1 <- fc[rn,]
      gsva1 <- gsva[rn,]

      name <- temp[1,"name"]
      x1 <- temp$bmd.fc
      x2 <- temp$bmd.gsva

      y1 <- temp$hitcall.fc
      y2 <- temp$hitcall.gsva

      plot(gsva1$top_over_cutoff~fc1$top_over_cutoff,xlab="SN(fc)",ylab=paste0("SN(",method2,")"),main=name,cex.lab=1.2,cex.axis=1.2,
           xlim=c(0,15),ylim=c(0,15))

      lines(c(0,100),c(0,100))
      if(!to.file) browser()
    }
  }

  fc$proper_name <- fc$name
  gsva$proper_name <- gsva$name

  mat <- mat[order(mat$delta,decreasing=T),]
  mat <- mat[mat$delta>100,]
  for(i in 1:nrow(mat)) {
    rn <- rownames(mat)[i]
    row.fc <- fc[rn,]
    row.gsva <- gsva[rn,]
    pathwayConcRespPlot(row.fc,CYTOTOX)
    pathwayConcRespPlot(row.gsva,CYTOTOX)
    if(!to.file) browser()
  }
  if(to.file) dev.off()
}

