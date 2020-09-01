library(openxlsx)
library(tidyr)
library(stringr)
library(stringi)
library(reshape2)
library(gplots)
library(RColorBrewer)
#--------------------------------------------------------------------------------------
#' Plot the results of refchem.replicability.stats
#'
#' heparg2d_toxcast_pfas_pe1_normal_refchems
#' u2os_toxcast_pfas_pe1_normal_refchems
#--------------------------------------------------------------------------------------
refchem.replicability.plots <- function(to.file=F,
                                        sigset="screen_large",
                                        method="fc") {
  printCurrentFunction(paste(sigset,method))

  dataset.1 <- "u2os_toxcast_pfas_pe1_normal_refchems"
  dataset.2 <- "heparg2d_toxcast_pfas_pe1_normal_refchems"
  file <- paste0("../output/signature_replicability/refchem.replicability.stats ",dataset.1,"_",sigset,"_",method,".xlsx")
  mat.1 <- read.xlsx(file)
  file <- paste0("../output/signature_replicability/refchem.replicability.stats ",dataset.2,"_",sigset,"_",method,".xlsx")
  mat.2 <- read.xlsx(file)
  mat <- rbind(mat.1,mat.2)
  if(to.file) {
    fname <- paste0("../output/signature_replicability/refchem.replicability.plots ",sigset,"_",method,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,1))
  ctype.list <- unique(mat$celltype)
  for(celltype in ctype.list) {
    matc <- mat[is.element(mat$celltype,celltype),]
    chem.list <- sort(unique(matc$chemical))
    for(chemical in chem.list) {
      matcc <- matc[is.element(matc$chemical,chemical),]

      ################################################################################
      ################################################################################
      ################################################################################
      ################################################################################
      x <- matcc$tc_mean
      y1 <- matcc$f_top_02
      y2 <- matcc$f_top_max_sign
      y3 <- matcc$f_bmd_10
      lims <- seq(from=1,to=4,by=0.25)

      ################################################################################
      cats <- NULL
      counts <- NULL
      for(lim in lims) {
        ytemp <- y1[x>lim]
        xtemp <- ytemp
        xtemp[] <- paste0(">",lim)
        cats <- c(cats,xtemp)
        counts <- c(counts,ytemp)
      }
      boxplot(counts~cats,main=paste("f(top) :",celltype,":",chemical),ylim=c(0,1),cex.lab=1.5,cex.axis=1.5,
              xlab="T/C",ylab="fraction")
      for(t in c(0.2,0.4,0.6,0.8,1.0)) lines(c(0,100),c(t,t))

      ################################################################################
      cats <- NULL
      counts <- NULL
      for(lim in lims) {
        ytemp <- y2[x>lim]
        xtemp <- ytemp
        xtemp[] <- paste0(">",lim)
        cats <- c(cats,xtemp)
        counts <- c(counts,ytemp)
      }
      boxplot(counts~cats,main=paste("f(sign) :",celltype,":",chemical),ylim=c(0,1),cex.lab=1.5,cex.axis=1.5,
              xlab="T/C",ylab="fraction")
      for(t in c(0.2,0.4,0.6,0.8,1.0)) lines(c(0,100),c(t,t))

      ################################################################################
      cats <- NULL
      counts <- NULL
      for(lim in lims) {
        ytemp <- y3[x>lim]
        xtemp <- ytemp
        xtemp[] <- paste0(">",lim)
        cats <- c(cats,xtemp)
        counts <- c(counts,ytemp)
      }
      boxplot(counts~cats,main=paste("f(bmd) :",celltype,":",chemical),ylim=c(0,1),cex.lab=1.5,cex.axis=1.5,
              xlab="T/C",ylab="fraction")
      for(t in c(0.2,0.4,0.6,0.8,1.0)) lines(c(0,100),c(t,t))

      if(!to.file) browser()

      ################################################################################
      ################################################################################
      ################################################################################
      ################################################################################
      x <- matcc$hc_mean
      y1 <- matcc$f_top_02
      y2 <- matcc$f_top_max_sign
      y3 <- matcc$f_bmd_10
      lims <- seq(from=0,to=1,by=0.1)

      ################################################################################
      cats <- NULL
      counts <- NULL
      for(lim in lims) {
        ytemp <- y1[x>lim]
        xtemp <- ytemp
        xtemp[] <- paste0(">",lim)
        cats <- c(cats,xtemp)
        counts <- c(counts,ytemp)
      }
      boxplot(counts~cats,main=paste("f(top) :",celltype,":",chemical),ylim=c(0,1),cex.lab=1.5,cex.axis=1.5,
              xlab="hitcall",ylab="fraction")
      for(t in c(0.2,0.4,0.6,0.8,1.0)) lines(c(0,100),c(t,t))

      ################################################################################
      cats <- NULL
      counts <- NULL
      for(lim in lims) {
        ytemp <- y2[x>lim]
        xtemp <- ytemp
        xtemp[] <- paste0(">",lim)
        cats <- c(cats,xtemp)
        counts <- c(counts,ytemp)
      }
      boxplot(counts~cats,main=paste("f(sign) :",celltype,":",chemical),ylim=c(0,1),cex.lab=1.5,cex.axis=1.5,
              xlab="hitcall",ylab="fraction")
      for(t in c(0.2,0.4,0.6,0.8,1.0)) lines(c(0,100),c(t,t))

      ################################################################################
      cats <- NULL
      counts <- NULL
      for(lim in lims) {
        ytemp <- y3[x>lim]
        xtemp <- ytemp
        xtemp[] <- paste0(">",lim)
        cats <- c(cats,xtemp)
        counts <- c(counts,ytemp)
      }
      boxplot(counts~cats,main=paste("f(bmd) :",celltype,":",chemical),ylim=c(0,1),cex.lab=1.5,cex.axis=1.5,
              xlab="hitcall",ylab="fraction")
      for(t in c(0.2,0.4,0.6,0.8,1.0)) lines(c(0,100),c(t,t))

      if(!to.file) browser()
    }
  }
  par(mfrow=c(6,2))
  ctype.list <- unique(mat$celltype)
  for(celltype in ctype.list) {
    matc <- mat[is.element(mat$celltype,celltype),]
    chem.list <- sort(unique(matc$chemical))
    for(chemical in chem.list) {
      matcc <- matc[is.element(matc$chemical,chemical),]
      breaks <- seq(from=0,to=1,by=0.1)
      hist(matcc$hc_mean,main=paste("Hitcall :",celltype,":",chemical),cex.lab=1.5,cex.axis=1.5,breaks=breaks,xlab="Hitcall")

      breaks <- seq(from=0,to=4,by=0.25)
      x <- matcc$tc_mean
      x <- x[x<4]
      hist(x,main=paste("T/C :",celltype,":",chemical),cex.lab=1.5,cex.axis=1.5,breaks=breaks,xlab="Top/Cutoff")

      if(!to.file) browser()
    }
  }

  if(to.file) dev.off()
}
