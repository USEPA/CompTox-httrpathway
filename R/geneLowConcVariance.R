library(openxlsx)
library(tidyr)
library(stringr)
library(stringi)
library(reshape2)
library(gplots)
library(RColorBrewer)
#--------------------------------------------------------------------------------------
#' Look for genes with high variance at low concentrations
#' @param method signature scoring method in c("fc", "gsva", "mygsea")
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#' u2os_meanncnt0_5-plateteffect_1-shrinkage_normal
#--------------------------------------------------------------------------------------
geneLowConcVariance <- function(to.file=F,
                                do.load=F,
                                dataset="mcf7_ph1_pe1_normal_good_pg",
                                celltype="MCF7") {
  printCurrentFunction(paste(dataset))
  if(to.file) {
    fname <- paste0("../input/fcdata/geneLowConcVariance ",dataset,".pdf")
    pdf(file=fname,width=6,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(2,1),mar=c(4,4,2,2))

  if(do.load) {
    cat("do.load\n")
    file <- paste0("../input/fcdata/FCMAT2_",dataset,".RData")
    print(file)
    load(file=file)
    FCMAT2 <<- FCMAT2
    file <- paste0("../input/fcdata/CHEM_DICT_",dataset,".RData")
    print(file)
    load(file=file)
    CHEM_DICT <<- CHEM_DICT
  }

  index.low <- CHEM_DICT[CHEM_DICT$conc_index=="1","sample_key"]
  index.high <- CHEM_DICT[CHEM_DICT$conc_index=="8","sample_key"]

  mat.low <- FCMAT2[index.low,]
  mat.high <- FCMAT2[index.high,]

  med.low <- apply(mat.low,FUN=median,MARGIN=2)
  mad.low <- apply(mat.low,FUN=mad,MARGIN=2)

  med.high <- apply(mat.high,FUN=median,MARGIN=2)
  mad.high <- apply(mat.high,FUN=mad,MARGIN=2)

  plot(med.low~med.high,cex.lab=1.2,cex.axis=1.2,
       main=paste0(celltype,": Median"),xlab="med(high conc)",ylab="med(low conc)",pch=".",
       xlim=c(-0.2,0.2),ylim=c(-0.2,0.2))
  lines(c(-1,2),c(-1,2))
  plot(mad.low~mad.high,cex.lab=1.2,cex.axis=1.2,
       main=paste0(celltype,": MAD"),xlab="mad(high conc)",ylab="mad(low conc)",pch=".",
       xlim=c(0,0.25),ylim=c(0,0.25))
  lines(c(0,2),c(0,2))

  res <- t(rbind(med.low,med.high,mad.low,mad.high))
  res <- res[order(res[,"mad.low"],decreasing=T),]
  res <- as.data.frame(res)
  res <- cbind(rownames(res),res)
  names(res)[1] <- "gene"

  file <- paste0("../input/fcdata/geneLowConcVariance ",dataset,".xlsx")
  write.xlsx(res,file)

  if(!to.file) browser()
  else dev.off()
}
