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
#' u2os_meanncnt0_5-plateteffect_1-shrinkage_normal
#--------------------------------------------------------------------------------------
dup.replicability.stats <- function(to.file=F,
                                    do.load=F,
                                    do.prep=F,
                                    dataset="mcf7_ph1_pe1_normal_good_pg",
                                    celltype="MCF7",
                                    sigset="screen_large",
                                    method="fc") {
  printCurrentFunction(paste(dataset,sigset,method))
  if(to.file) {
    fname <- paste0("../output/signature_replicability/dup.replicability.stats ",dataset,"_",sigset,"_",method,".pdf")
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

  if(do.prep) {
    cat("  do.prep 1\n")
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
    cat("  do.prep 2\n")

    file <- "../input/pfas/QC sample map 2020-05-04.xlsx"
    qc <- read.xlsx(file)
    qc <- qc[!is.element(qc$score,c("H","M")),]
    spid.list <- qc$spid
    mat <- mat[!is.element(mat$sample_id,spid.list),]
    cat("  do.prep 3\n")

    er <- mat$er
    cutoff <- mat$cutoff
    err.range <- exp(er)*2.7765
    err_over_cutoff <- err.range / mat$cutoff
    mat$err_over_cutoff <- err_over_cutoff
    cat("  do.prep 4\n")

    dtxsid.list <- unique(mat$dtxsid)
    chems <- unique(mat[,c("dtxsid","sample_id")])
    chems <- chems[order(chems$dtxsid),]
    chems$set <- 0
    chems[1,"set"] <- 1
    for(i in 2:nrow(chems)) {
      if(chems[i,"dtxsid"]==chems[i-1,"dtxsid"]) {
        if(chems[i-1,"set"]==1) chems[i,"set"] <- 2
        else if(chems[i-1,"set"]==2) chems[i,"set"] <- 3
      }
      else chems[i,"set"] <- 1
    }
    cat("  do.prep 5\n")

    chems <- chems[chems$set>0,]
    sid1 <- chems[chems$set==1,"sample_id"]
    sid2 <- chems[chems$set==2,"sample_id"]
    mat1 <- mat[is.element(mat$sample_id,sid1),]
    mat2 <- mat[is.element(mat$sample_id,sid2),]

    mat1 <- mat1[order(mat1$signature),]
    mat2 <- mat2[order(mat2$signature),]
    cat("  do.prep 6\n")

    res <- cbind(mat1[c("dtxsid","casrn","name","signature",
                        "sample_id","top_over_cutoff","hitcall","err_over_cutoff","bmd","top")],
                 mat2[c("sample_id","top_over_cutoff","hitcall","err_over_cutoff","bmd","top")])
    name.list <- c("dtxsid","casrn","name","signature",
                   "sid1","tc1","hc1","ec1","bmd1","top1",
                   "sid2","tc2","hc2","ec2","bmd2","top2")

    names(res) <- name.list

    cat("  do.prep 7\n")
    res$delta_bmd <- abs(log10(res$bmd1)-log10(res$bmd2))
    delta_tsign <- res$top1/abs(res$top1) - res$top2/abs(res$top2)
    delta_tsign <- -abs(delta_tsign)/2
    delta_tsign[delta_tsign==0] <- 1
    res$delta_top <- delta_tsign * abs(res$top1-res$top2)

    res[is.na(res$tc1),"tc1"] <- 0
    res[is.na(res$tc2),"tc2"] <- 0
    res[is.na(res$hc1),"hc1"] <- 0
    res[is.na(res$hc2),"hc2"] <- 0
    res[is.na(res$ec1),"ec1"] <- 0
    res[is.na(res$ec2),"ec2"] <- 0

    res$max_tc <- res$tc1
    res[res$tc2>res$tc1,"max_tc"] <- res[res$tc2>res$tc1,"tc2"]

    res$max_hc <- res$hc1
    res[res$hc2>res$hc1,"max_hc"] <- res[res$hc2>res$hc1,"hc2"]

    res$min_ec <- res$ec1
    res[res$ec2<res$ec1,"max_ec"] <- res[res$ec2<res$ec1,"ec2"]

    RES <<- res
    cat("  do.prep 8\n")

    file <- paste0("../output/signature_replicability/dup.replicability.stats ",dataset,"_",sigset,"_",method,".RData")
    save(RES,file=file)
    cat("  do.prep 9\n")
  }
  res <- RES
  ####################################################################
  ####################################################################
  ####################################################################
  ####################################################################
  ####################################################################
  ####################################################################
  ####################################################################

  hc.list <- c(0.5,0.6,0.7,0.8,0.9)
  x <- NULL
  y <- NULL
  dbmd.list <- hc.list
  dbmd.list[] <- 0
  for(i in 1:length(hc.list)) {
    hc <- hc.list[i]
    dbmd <- res[res$max_hc>hc,"delta_bmd"]
    temp <- vector(length=length(dbmd),mode="character")
    temp[] <- paste0(">",hc)
    x <- c(x,temp)
    y <- c(y,dbmd)
    dbmd.list[i] <- length(dbmd[dbmd<1])/length(dbmd)
  }
  plot(dbmd.list~hc.list,type="h",lwd=10,xlim=c(0.5,0.9),ylim=c(0,1),xlab="Hitcall > x",ylab="f(|BMD1-BMD2|<10)",cex.lab=1.5,cex.axis=1.5,main=celltype)
  for(top in c(0.2,0.4,0.6,0.8,1.0)) lines(c(0,100),c(top,top))
  boxplot(y~x,main="delta(bmd)~HC",ylim=c(0,4),xlab="HC Cutoff",ylab="delta log(BMD)")

  x <- NULL
  y <- NULL
  dtop.list <- hc.list
  dtop.list[] <- 0
  for(i in 1:length(hc.list)) {
    hc <- hc.list[i]
    dtop <- res[res$max_hc>hc,"delta_top"]
    temp <- vector(length=length(dtop),mode="character")
    temp[] <- paste0(">",hc)
    x <- c(x,temp)
    y <- c(y,dtop)
    denom <- length(dtop)
    dtop <- dtop[dtop>0]
    dtop <- dtop[dtop<0.2]
    dtop.list[i] <- length(dtop)/denom
  }
  plot(dtop.list~hc.list,type="h",lwd=10,xlim=c(0.5,0.9),ylim=c(0,1),xlab="Hitcall > x",ylab="f(|Top1-Top2|<0.2)",cex.lab=1.5,cex.axis=1.5,main=celltype)
  for(top in c(0.2,0.4,0.6,0.8,1.0)) lines(c(0,100),c(top,top))
  boxplot(y~x,main="delta(top)~HC",ylim=c(-1,1),xlab="HC Cutoff",ylab="delta log(Top)")
  #####################################################################################

  res <- res[res$max_hc>0.8,]
  tc.list <- seq(from=1.0,to=4.0,by=0.5)
  x <- NULL
  y <- NULL
  dbmd.list <- tc.list
  dbmd.list[] <- 0
  for(i in 1:length(tc.list)) {
    tc <- tc.list[i]
    dbmd <- res[res$max_tc>tc,"delta_bmd"]
    temp <- vector(length=length(dbmd),mode="character")
    temp[] <- paste0(">",tc)
    x <- c(x,temp)
    y <- c(y,dbmd)
    dbmd.list[i] <- length(dbmd[dbmd<1])/length(dbmd)
  }
  plot(dbmd.list~tc.list,type="h",lwd=10,xlim=c(1,4),ylim=c(0,1),xlab="Top/Cutoff > x",ylab="f(|BMD1-BMD2|<10)",cex.lab=1.5,cex.axis=1.5,main=celltype)
  for(top in c(0.2,0.4,0.6,0.8,1.0)) lines(c(0,100),c(top,top))
  boxplot(y~x,main="delta(bmd)~TC",ylim=c(0,4),xlab="TC Cutoff",ylab="delta log(BMD)")

  x <- NULL
  y <- NULL
  dtop.list <- tc.list
  dtop.list[] <- 0
  for(i in 1:length(tc.list)) {
    tc <- tc.list[i]
    dtop <- res[res$max_tc>tc,"delta_top"]
    temp <- vector(length=length(dtop),mode="character")
    temp[] <- paste0(">",tc)
    x <- c(x,temp)
    y <- c(y,dtop)
    denom <- length(dtop)
    dtop <- dtop[dtop>0]
    dtop <- dtop[dtop<0.2]
    dtop.list[i] <- length(dtop)/denom
    #browser()
  }
  plot(dtop.list~tc.list,type="h",lwd=10,xlim=c(1,4),ylim=c(0,1),xlab="Top/Cutoff > x",ylab="f(|Top1-Top2|<0.2)",cex.lab=1.5,cex.axis=1.5,main=celltype)
  for(top in c(0.2,0.4,0.6,0.8,1.0)) lines(c(0,100),c(top,top))
  boxplot(y~x,main="delta(top)~TC",ylim=c(-1,1),xlab="TC Cutoff",ylab="delta log(Top)")
  if(!to.file) browser()

  ################################################################################
  res <- res[res$max_tc>2,]
  ec.list <- seq(from=0,to=10,by=2)
  x <- NULL
  y <- NULL
  for(ec in ec.list) {
    dbmd <- res[res$max_ec<ec,"delta_bmd"]
    temp <- vector(length=length(dbmd),mode="character")
    temp[] <- paste0("<",ec)
    x <- c(x,temp)
    y <- c(y,dbmd)
  }
  boxplot(y~x,main="delta(bmd)~EC",ylim=c(0,4),xlab="EC Cutoff",ylab="delta log(BMD)")

  x <- NULL
  y <- NULL
  for(ec in ec.list) {
    dbmd <- res[res$max_ec<ec,"delta_top"]
    temp <- vector(length=length(dbmd),mode="character")
    temp[] <- paste0("<",ec)
    x <- c(x,temp)
    y <- c(y,dbmd)
  }
  boxplot(y~x,main="delta(top)~EC",ylim=c(-1,1),xlab="EC Cutoff",ylab="delta log(Top)")

  if(!to.file) browser()
  ####################################################################
  ####################################################################
  ####################################################################
  ####################################################################
  ####################################################################
  ####################################################################
  ####################################################################
  ####################################################################


  hc.list <- c(0.5,0.6,0.7,0.8,0.9)
  x <- NULL
  y <- NULL
  for(hc in hc.list) {
    dbmd <- res[res$max_hc>hc,"delta_bmd"]
    temp <- vector(length=length(dbmd),mode="character")
    temp[] <- paste0(">",hc)
    x <- c(x,temp)
    y <- c(y,dbmd)
  }
  boxplot(y~x,main="delta(bmd)~HC",ylim=c(0,4),xlab="HC Cutoff",ylab="delta log(BMD)")

  x <- NULL
  y <- NULL
  for(hc in hc.list) {
    dbmd <- res[res$max_hc>hc,"delta_top"]
    temp <- vector(length=length(dbmd),mode="character")
    temp[] <- paste0(">",hc)
    x <- c(x,temp)
    y <- c(y,dbmd)
  }
  boxplot(y~x,main="delta(top)~HC",ylim=c(-1,1),xlab="HC Cutoff",ylab="delta log(Top)")
  #####################################################################################

  res <- res[res$max_hc>0.8,]
  tc.list <- seq(from=1.0,to=4.0,by=0.5)
  x <- NULL
  y <- NULL
  for(tc in tc.list) {
    dbmd <- res[res$max_tc>tc,"delta_bmd"]
    temp <- vector(length=length(dbmd),mode="character")
    temp[] <- paste0(">",tc)
    x <- c(x,temp)
    y <- c(y,dbmd)
  }
  boxplot(y~x,main="delta(bmd)~TC",ylim=c(0,4),xlab="TC Cutoff",ylab="delta log(BMD)")

  x <- NULL
  y <- NULL
  for(tc in tc.list) {
    dbmd <- res[res$max_tc>tc,"delta_top"]
    temp <- vector(length=length(dbmd),mode="character")
    temp[] <- paste0(">",tc)
    x <- c(x,temp)
    y <- c(y,dbmd)
  }
  boxplot(y~x,main="delta(top)~TC",ylim=c(-1,1),xlab="TC Cutoff",ylab="delta log(Top)")

  ################################################################################
  res <- res[res$max_tc>2,]
  ec.list <- seq(from=0,to=10,by=2)
  x <- NULL
  y <- NULL
  for(ec in ec.list) {
    dbmd <- res[res$max_ec<ec,"delta_bmd"]
    temp <- vector(length=length(dbmd),mode="character")
    temp[] <- paste0("<",ec)
    x <- c(x,temp)
    y <- c(y,dbmd)
  }
  boxplot(y~x,main="delta(bmd)~EC",ylim=c(0,4),xlab="EC Cutoff",ylab="delta log(BMD)")

  x <- NULL
  y <- NULL
  for(ec in ec.list) {
    dbmd <- res[res$max_ec<ec,"delta_top"]
    temp <- vector(length=length(dbmd),mode="character")
    temp[] <- paste0("<",ec)
    x <- c(x,temp)
    y <- c(y,dbmd)
  }
  boxplot(y~x,main="delta(top)~EC",ylim=c(-1,1),xlab="EC Cutoff",ylab="delta log(Top)")

  if(!to.file) browser()
  ####################################################################
  ####################################################################
  ####################################################################
  ####################################################################
  ####################################################################
  ####################################################################
  ####################################################################
  ####################################################################
  cat("  do cr plots 1\n")

  par(mfrow=c(3,2),mar=c(4,4,2,2))
  res <- RES

  res <- res[order(res$max_tc,decreasing=T),]
  mask <- vector(mode="integer",length=nrow(res))
  mask[] <- 0
  mask1 <- mask;mask2 <- mask
  mask1[res$tc1>2] <- 1
  mask2[res$tc2>2] <- 1
  mask12 <- mask1+mask2
  res <- res[mask12>0,]
  cat("  do cr plots 2a:",nrow(res),"\n")

  mask <- vector(mode="integer",length=nrow(res))
  mask1 <- mask;mask2 <- mask
  mask1[res$hc1>0.7] <- 1
  mask2[res$hc2>0.7] <- 1
  mask12 <- mask1+mask2
  res <- res[mask12>0,]
  cat("  do cr plots 2b:",nrow(res),"\n")

  mask <- vector(mode="integer",length=nrow(res))
  mask1 <- mask;mask2 <- mask
  mask1[res$bmd1<10] <- 1
  mask2[res$bmd2<10] <- 1
  mask12 <- mask1+mask2
  res <- res[mask12>0,]
  cat("  do cr plots 3:",nrow(res),"\n")

  list1 <- paste0(res$sid1,"_",res$signature)
  list2 <- paste0(res$sid2,"_",res$signature)
  d1 <- as.data.frame(matrix(nrow=length(list1),ncol=2))
  d2 <- d1
  d1[,1] <- list1
  d2[,1] <- list2
  n <- length(list1)
  vals1 <- seq(from=1,to=n*2,by=2)
  vals1 <- vals1[1:n]
  vals2 <- vals1+1
  d1[,2] <- vals1
  d2[,2] <- vals2
  cat("  do cr plots 4\n")

  dmat <- rbind(d1,d2)
  dmat <- dmat[order(dmat[2]),]
  row.list <- dmat[,1]
  cat("  do cr plots 5:",length(row.list),"\n")

  mat <- MAT
  sid.list <- unique(c(res$sid1,res$sid2))
  mat <- mat[is.element(mat$sample_id,sid.list),]
  mat <- mat[is.element(mat$signature,res$signature),]
  cat("  do cr plots 6\n")

  rnames <- paste0(mat$sample_id,"_",mat$signature)
  cat("  do cr plots 7\n")

  rownames(mat) <- rnames
  cat("  do cr plots 8\n")

  mat1 <- mat[row.list,]
  cat(" number of rows:",nrow(mat1),"\n")

  mat1$proper_name <- mat1$name
  cat("  cycle through signatures (rows) and run signatureConcRespPlot\n")

  for(i in 1:nrow(mat1)){
    cat(i," out of ",nrow(mat1),"\n")
    signatureConcRespPlot(mat1[i,])
    if(!to.file) browser()
  }

  if(!to.file) browser()
  else dev.off()
}
