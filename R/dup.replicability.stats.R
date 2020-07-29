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
dup.replicability.stats <- function(to.file=F,
                                    do.prep=F,
                                    do.load=F,
                                    dataset="heparg2d_toxcast_pfas_pe1_normal",
                                    sigcatalog="signatureDB_master_catalog 2020-07-10",
                                    sigset="screen_large",
                                    method="mygsea") {
  printCurrentFunction(paste(dataset,sigset,method))
  if(to.file) {
    fname <- paste0("../output/signature_replicability/dup.replicability.stats ",dataset,"_",sigset,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,2,2))

  if(do.load) {
    cat("do.load\n")
    annotations <- signatureCatalogLoader(sigset,sigcatalog)
    ANN <<- annotations

    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat <- SIGNATURE_CR
    MAT <<- mat
  }

  if(do.prep) {
    cat("do.prep\n")
    mat <- MAT
    mat[is.na(mat$super_target),"super_target"] <- "-"
    mat <- mat[!is.element(mat$super_target,"-"),]
    catalog <- ANN
    res1 <- unique(mat[,c("sample_id","dtxsid","casrn","name")])
    res1 <- res1[order(res1$dtxsid),]
    rownames(res1) <- res1$sample_id
    dtxsid.list <- res1$dtxsid
    dups <- duplicated(dtxsid.list)
    temp <- dtxsid.list[dups]
    res1$duplicate <- 0
    res1[is.element(res1$dtxsid,temp),"duplicate"] <- 1
    dups <- res1[res1$duplicate==1,]
    mat <- mat[is.element(mat$sample_id,dups$sample_id),]

    file <- "../input/pfas/QC sample map 2020-05-04.xlsx"
    qc <- read.xlsx(file)
    qc <- qc[!is.element(qc$score,c("H","M")),]
    spid.list <- qc$spid
    mat <- mat[!is.element(mat$sample_id,spid.list),]

    er <- mat$er
    cutoff <- mat$cutoff
    err.range <- exp(er)*2.7765
    err_over_cutoff <- err.range / mat$cutoff
    mat$err_over_cutoff <- err_over_cutoff

    dtxsid.list <- unique(mat$dtxsid)
    chems <- unique(mat[,c("dtxsid","sample_id")])
    chems <- chems[order(chems$dtxsid),]
    chems$set <- 0
    chems[1,"set"] <- 1
    for(i in 2:nrow(chems)) {
      if(chems[i,"dtxsid"]==chems[i-1,"dtxsid"]) {
        if(chems[i-1,"set"]==1) chems[i,"set"] <- 2
      }
      else chems[i,"set"] <- 1
    }
    chems <- chems[chems$set>0,]
    sid1 <- chems[chems$set==1,"sample_id"]
    sid2 <- chems[chems$set==2,"sample_id"]
    mat1 <- mat[is.element(mat$sample_id,sid1),]
    mat2 <- mat[is.element(mat$sample_id,sid2),]

    mat1 <- mat1[order(mat1$signature),]
    mat2 <- mat2[order(mat2$signature),]

    res <- cbind(mat1[c("dtxsid","casrn","name","signature",
                        "sample_id","top_over_cutoff","hitcall","err_over_cutoff","bmd")],
                 mat2[c("sample_id","top_over_cutoff","hitcall","err_over_cutoff","bmd")])
    name.list <- c("dtxsid","casrn","name","signature",
                   "sid1","tc1","hc1","ec1","bmd1",
                   "sid2","tc2","hc2","ec2","bmd2")

    names(res) <- name.list
    RES <<- res

    file <- paste0("../output/signature_replicability/dup.replicability.stats ",dataset,"_",sigset,".RData")
    save(RES,file=file)
  }
  res <- RES

  hc.list <- c(0.9,0.8,0.7,0.6,0.5)
  tc.list <- c(5,4.5,4,3.5,3,2.5,2,1.5,1)
  ec.list <- c(0)
  name.list <- c("hc","tc","ec","frac")
  perc <- matrix(nrow=length(hc.list)*length(tc.list)*length(ec.list),ncol=length(name.list))
  colnames(perc) <- name.list
  counter <- 0
  for(ec in ec.list) {
    for(hc in hc.list) {
      for(tc in tc.list) {
        counter <- counter + 1
        ecmin <- ec
        ecmax <- ecmin+6
        tcmin <- tc
        tcmax <- tcmin+0.5
        hcmin <- hc
        hcmax <- hcmin+0.1
        temp <- res[res$hc1>hcmin,]
        temp <- temp[temp$hc1<=hcmax,]

        temp <- temp[temp$tc1>tcmin,]
        temp <- temp[temp$tc1<=tcmax,]

        temp <- temp[temp$ec1>ecmin,]
        temp <- temp[temp$ec1<=ecmax,]

        perc[counter,"hc"] <- hc
        perc[counter,"tc"] <- tc
        perc[counter,"ec"] <- ec
        if(nrow(temp)<=3) perc[counter,"frac"] <- 0
        else {
          ngt5 <- nrow(temp[temp$hc2>0.5,])
          perc[counter,"frac"] <- ngt5 / nrow(temp)
        }
      }
    }
  }
  for(ec in ec.list) {
    plot(c(0,0),xlab="hc",ylab="tc",xlim=c(0.5,1),ylim=c(1,5),cex.lab=1.2,cex.axis=1.2,type="n",
         main=paste("E/C:",ec))
    temp <- perc[perc[,"ec"]==ec,]
    for(i in 1:nrow(temp)) {
      x <- temp[i,"hc"]
      y <- temp[i,"tc"]
      f <- temp[i,"frac"]
      color <- "black"
      if(f>0.5) color <- "red"
      if(f>0) text(x,y,format(f,digits=2),cex=0.9,col=color)
    }
    if(!to.file) browser()
  }

  plot(res$tc2~res$tc1,cex.lab=1.2,cex.axis=1.2,main="TC",xlim=c(0,6),ylim=c(0,6),pch=".")
  plot(res$hc2~res$hc1,cex.lab=1.2,cex.axis=1.2,main="HC",xlim=c(0,1),ylim=c(0,1),pch=".")
  plot(res$ec2~res$ec1,cex.lab=1.2,cex.axis=1.2,main="EC",xlim=c(0,10),ylim=c(0,10),pch=".")

  par(mfrow=c(3,2),mar=c(4,4,2,2))
  res1 <- res[res$hc2>0.5,]
  res2 <- res[res$hc2<=0.5,]
  hist(res1$hc1)
  hist(res2$hc1)

  hist(res1$tc1)
  hist(res2$tc1)

  hist(res1$ec1)
  hist(res2$ec1)

  plot(density(res1$hc1,na.rm=T))
  lines(density(res2$hc1,na.rm=T),col="red")

  plot(density(res1$tc1,na.rm=T))
  lines(density(res2$tc1,na.rm=T),col="red")

  plot(density(res1$ec1,na.rm=T))
  lines(density(res2$ec1,na.rm=T),col="red")

  plot(res1$hc1~res1$ec1,pch=".",col="black")
  points(res2$hc1~res2$ec1,pch=".",col="red")
  dtxsid.list <- unique(res1$dtxsid)
  dtxsid.list <- dtxsid.list[is.element(dtxsid.list,res2$dtxsid)]
  for(dtxsid in dtxsid.list) {
    res1d <- res1[is.element(res1$dtxsid,dtxsid),]
    res2d <- res2[is.element(res2$dtxsid,dtxsid),]
    name <- res1d[1,"name"]
    plot(res1d$hc1~res1d$tc1,pch=".",col="black",main=name,xlim=c(0,6),ylim=c(0,1))
    points(res2d$hc1~res2d$tc1,pch=".",col="red")
  }
  if(!to.file) browser()
  else dev.off()
}
