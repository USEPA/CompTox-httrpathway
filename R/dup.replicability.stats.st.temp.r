library(openxlsx)
library(tidyr)
library(stringr)
library(stringi)
library(reshape2)
library(gplots)
library(RColorBrewer)
#--------------------------------------------------------------------------------------
#' Compare pairs of duplicated chemical / signatures to see what parameter
#' drive replicability
#' @param method signature scoring method in c("fc", "gsva", "mygsea")
#'
#'
#' Error bars are exp(er)*qt(.025,4) = exp(er)*2.7765
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#--------------------------------------------------------------------------------------
dup.replicability.stats.st.temp <- function(to.file=F,
                                            do.load=F,
                                            dataset="heparg2d_toxcast_pfas_pe1_normal",
                                            sigset="screen_large") {
  printCurrentFunction(paste(dataset,sigset))
  if(to.file) {
    fname <- paste0("../output/signature_replicability/dup.replicability.stats.st.temp ",dataset,"_",sigset,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,2,2))

  if(do.load) {
    cat("do.load\n")
    file <- paste0("../output/signature_replicability/dup.replicability.stats.st res2 ",dataset,"_",sigset,".RData")
    print(file)
    load(file=file)

    RES2 <<- RES2
  }
  res2 <- RES2
  plot(res2$tc2~res2$tc1,pch=".")
  plot(res2$hc2~res2$hc1,pch=".")
  plot(res2$ec2~res2$ec1,pch=".")
  plot(res2$bmd2~res2$bmd1,log="xy",pch=".")
  plot(res2$caikwt2~res2$caikwt1,pch=".")

  if(!to.file) browser()
  else dev.off()
}
