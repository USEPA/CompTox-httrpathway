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
#' heparg2d_toxcast_pfas_pe1_normal_refchems
#' u2os_toxcast_pfas_pe1_normal_refchems
#--------------------------------------------------------------------------------------
refchem.replicability.stats <- function(to.file=F,
                                        do.load=F,
                                        dataset="heparg2d_toxcast_pfas_pe1_normal_refchems",
                                        celltype="HepaRG",
                                        method="fc",
                                        sigset="screen_large") {
  printCurrentFunction(paste(dataset,sigset,method))

  if(do.load) {
    cat("do.load\n")
    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat <- SIGNATURE_CR

    MAT <<- mat
  }
  mat <- MAT

  t1 <- mat$sample_id
  t3 <- str_replace(t1,"BAP.","")
  t3 <- str_replace(t3,"RIF","")
  t3 <- str_replace(t3,"TROG.","")
  mat$pg <- t3
  if(celltype=="HepaRG") {
    bad.set <- c("7","10","11","12","13","14")
    goopy.set <- c("3","8","9")
    cat(nrow(mat),"\n")
    mat <- mat[!is.element(mat$pg,bad.set),]
    cat(nrow(mat),"\n")
  }
  if(celltype=="U2OS") {
    bad.set <- c("19")
    browser()
  }
  file <- "../input/signatures/reference chemical super targets.xlsx"
  rcst <- read.xlsx(file)
  rcst <- rcst[!is.na(rcst$chemical),]
  rcst <- rcst[is.element(rcst$celltype,celltype),]

  mat <- mat[is.element(mat$super_target,rcst$super_target),]
  cat(nrow(mat),"\n")

  name.list <- c("celltype","chemical","signature","super_target","top_mean","top_sd","logbmd_mean","log_bmd_sd",
                 "tc_mean","tc_sd","hc_mean","hc_sd","f_top_max_sign","f_top_02","f_bmd_10","nrep")
  res.all <- NULL
  chem.list <- unique(mat$name)

  for(chem in chem.list) {
    if(to.file) {
      fname <- paste0("../output/signature_replicability/refchem.replicability.stats ",dataset,"_",sigset,"_",method,"_",chem,".pdf")
      pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    }
    par(mfrow=c(4,2),mar=c(4,4,2,2))

    matc <- mat[is.element(mat$name,chem),]

    rcstc <- rcst[is.element(rcst$chemical,chem),]
    cat("start:",chem,":",nrow(matc),"\n")
    matc <- matc[is.element(matc$super_target,rcstc$super_target),]
    cat("start:",chem,":",nrow(matc),"\n")
    temp <- matc[matc$top_over_cutoff>1,]
    temp <- temp[temp$hitcall>0.5,]
    #temp <- temp[temp$bmd<10,]
    slist <- sort(unique(temp$signature))
    nsig <- length(slist)
    res <- as.data.frame(matrix(nrow=nsig,ncol=length(name.list)))
    names(res) <- name.list
    res$celltype <- celltype
    res$chemical <- chem
    for(i in 1:nsig) {
      if(i%%100==0) cat(chem," ",i," out of ",nsig,"\n")
      signature <- slist[i]
      res[i,"signature"] <- signature
      temp <- matc[is.element(matc$signature,signature),]

      x <- temp$bmd
      y <- temp$top
      mask <- x
      mask[] <- 1
      mask[is.na(x)] <- 0
      mask[is.na(y)] <- 0
      x <- x[mask==1]
      y <- y[mask==1]

      res[i,"super_target"] <- temp[1,"super_target"]
      res[i,"top_mean"] <- mean(temp$top)
      res[i,"top_sd"] <- sd(temp$top)
      res[i,"logbmd_mean"] <- mean(log10(temp$bmd),na.rm=T)
      res[i,"log_bmd_sd"] <- sd(log10(temp$bmd),na.rm=T)
      res[i,"tc_mean"] <- mean(temp$top_over_cutoff,na.rm=T)
      res[i,"tc_sd"] <- sd(temp$top_over_cutoff,na.rm=T)
      res[i,"hc_mean"] <- mean(temp$hitcall,na.rm=T)
      res[i,"hc_sd"] <- sd(temp$hitcall,na.rm=T)
      res[i,"nrep"] <- length(x)

      fpos <- length(y[y>0])/length(y)
      fneg <- length(y[y<=0])/length(y)
      res[i,"f_top_max_sign"] <- max(fpos,fneg)
      mid <- mean(y)
      tin <- y[y>mid-0.1]
      tin <- tin[tin<mid+0.1]
      ymin <- mid-0.1
      ymax <- mid+0.1
      res[i,"f_top_02"] <- length(tin)/length(y)
      lx <- log10(x)
      mid <- mean(lx)
      bin <- lx[lx>mid-0.5]
      bin <- bin[bin<mid+0.5]
      xmin <- 10**(mid-0.5)
      xmax <- 10**(mid+0.5)
      res[i,"f_bmd_10"] <- length(bin)/length(x)

      plot(y~x,pch=".",xlab="BMD",ylab="Top",cex.lab=1.2,cex.axis=1.2,log="x",
           xlim=c(0.01,1000),ylim=c(-1,1),main=paste0(chem,"\n",signature),cex.main=0.9)
      lines(c(1e-5,1e5),c(0,0))
      text(1e-2,-1,res[i,"super_target"],pos=4)
      for(j in 1:length(x)) {
        xj <- x[j]
        yj <- y[j]
        col <- "black"
        if(xj<xmin || xj>xmax || yj<ymin || yj>ymax) col <- "red"
        points(xj,yj,pch=21,cex=1,bg=col,fg=col)
      }
      rect(xmin,ymin,xmax,ymax)
      if(!to.file) browser()
    }
    res.all <- rbind(res.all,res)
    if(to.file) dev.off()
  }
  file <- paste0("../output/signature_replicability/refchem.replicability.stats ",dataset,"_",sigset,"_",method,".xlsx")
  write.xlsx(res.all,file)
}
