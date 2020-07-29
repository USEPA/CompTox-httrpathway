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
dup.replicability.stats.st <- function(to.file=F,
                                       do.load=F,
                                       do.res1=F,
                                       do.res2=F,
                                       dataset="heparg2d_toxcast_pfas_pe1_normal",
                                       sigset="screen_large",
                                       method="mygsea") {
  printCurrentFunction(paste(dataset,sigset,method))
  if(to.file) {
    fname <- paste0("../output/signature_replicability/dup.replicability.stats.st ",dataset,"_",sigset,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,2,2))

  if(do.load) {
    cat("do.load\n")
    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat <- SIGNATURE_CR
    MAT <<- mat
  }

  if(do.res1) {
    cat("  do.res1\n")
    mat <- MAT
    mat[is.na(mat$super_target),"super_target"] <- "-"
    mat <- mat[!is.element(mat$super_target,"-"),]
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
        if(chems[i-1,"set"]==2) chems[i,"set"] <- 3
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

    res <- cbind(mat1[c("dtxsid","casrn","name","signature","super_target",
                        "sample_id","top_over_cutoff","hitcall","err_over_cutoff","bmd","caikwt")],
                 mat2[c("sample_id","top_over_cutoff","hitcall","err_over_cutoff","bmd","caikwt")])
    name.list <- c("dtxsid","casrn","name","signature","super_target",
                   "sid1","tc1","hc1","ec1","bmd1","caikwt1",
                   "sid2","tc2","hc2","ec2","bmd2","caikwt2")

    names(res) <- name.list
    RES1 <<- res

    file <- paste0("../output/signature_replicability/dup.replicability.stats.st res1 ",dataset,"_",sigset,".RData")
    save(RES1,file=file)
  }
  res <- RES1
  if(do.res2) {
    cat("  organize by super target\n")
    dtxsid.list <- unique(res$dtxsid)

    x <- unique(res[,c("signature","super_target")])
    st.list <- unique(x[duplicated(x$super_target),"super_target"])

    name.list <- names(res)
    ndt <- length(dtxsid.list)
    nst <- length(st.list)
    res2 <- as.data.frame(matrix(nrow=ndt*nst,ncol=length(name.list)))
    names(res2) <- name.list
    counter <- 0
    for(dtxsid in dtxsid.list) {
      cat(dtxsid,"\n")
      for(st in st.list) {
        temp <- res[is.element(res$dtxsid,dtxsid),]
        temp <- temp[is.element(temp$super_target,st),]
        temp <- temp[order(temp$tc1,decreasing=T),]
        res2[counter,1:11] <- temp[1,1:11]
        temp <- temp[order(temp$tc2,decreasing=T),]
        res2[counter,12:17] <- temp[1,12:17]
        counter <- counter+1
      }
    }
    res2$hit <- 0
    res2[is.na(res2$tc1),"tc1"] <- 0
    res2[is.na(res2$tc2),"tc2"] <- 0
    res2[is.na(res2$bmd1),"bmd1"] <- 1000
    res2[is.na(res2$bmd2),"bmd2"] <- 1000
    res2[res2$tc1>1,"hit"] <- 1
    res2[res2$tc2>1,"hit"] <- 1
    RES2 <<- res2
    file <- paste0("../output/signature_replicability/dup.replicability.stats.st res2 ",dataset,"_",sigset,".RData")
    save(RES2,file=file)
  }
  res2 <- RES2[RES2$hit==1,]

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
        temp <- res[res2$hc1>hcmin,]
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

  plot(res2$tc2~res2$tc1,cex.lab=1.2,cex.axis=1.2,main="TC",xlim=c(0,6),ylim=c(0,6),pch=".")
  plot(res2$hc2~res2$hc1,cex.lab=1.2,cex.axis=1.2,main="HC",xlim=c(0,1),ylim=c(0,1),pch=".")
  plot(res2$ec2~res2$ec1,cex.lab=1.2,cex.axis=1.2,main="EC",xlim=c(0,10),ylim=c(0,10),pch=".")

  par(mfrow=c(3,2),mar=c(4,4,2,2))
  resa <- res2[res$hc2>0.5,]
  resb <- res2[res$hc2<=0.5,]
  hist(resa$hc1)
  hist(resb$hc1)

  hist(resa$tc1)
  hist(resb$tc1)

  hist(resa$ec1)
  hist(resb$ec1)

  plot(density(resa$hc1,na.rm=T))
  lines(density(resb$hc1,na.rm=T),col="red")

  plot(density(resa$tc1,na.rm=T))
  lines(density(resb$tc1,na.rm=T),col="red")

  plot(density(resa$ec1,na.rm=T))
  lines(density(resb$ec1,na.rm=T),col="red")

  plot(resa$hc1~resa$ec1,pch=".",col="black")
  points(resb$hc1~resb$ec1,pch=".",col="red")
  dtxsid.list <- unique(resa$dtxsid)
  dtxsid.list <- dtxsid.list[is.element(dtxsid.list,resb$dtxsid)]
  for(dtxsid in dtxsid.list) {
    res1d <- resa[is.element(resa$dtxsid,dtxsid),]
    res2d <- resb[is.element(resb$dtxsid,dtxsid),]
    name <- res1d[1,"name"]
    plot(res1d$hc1~res1d$tc1,pch=".",col="black",main=name,xlim=c(0,6),ylim=c(0,1))
    points(res2d$hc1~res2d$tc1,pch=".",col="red")
  }
  if(!to.file) browser()
  else dev.off()
}
