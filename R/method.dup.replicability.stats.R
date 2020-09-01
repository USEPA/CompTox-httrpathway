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
method.dup.replicability.stats <- function(to.file=F,
                                           do.load=F,
                                           do.prep=F,
                                           dataset="mcf7_ph1_pe1_normal_good_pg",
                                           sigset="screen_large") {
  printCurrentFunction(paste(dataset,sigset))
  if(to.file) {
    fname <- paste0("../output/signature_replicability/method.dup.replicability.stats ",dataset,"_",sigset,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(6,4,2,2))

  if(do.load) {
    cat("   do.load\n")

    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_mygsea_0.05_conthits.RData")
    print(file)
    load(file=file)
    scr1 <- SIGNATURE_CR
    scr1[is.na(scr1$super_target),"super_target"] <- "-"
    scr1 <- scr1[!is.element(scr1$super_target,"-"),]
    rn1 <- paste0(scr1$sample_id,"_",scr1$signature)
    rownames(scr1) <- rn1

    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_fc_0.05_conthits.RData")
    print(file)
    load(file=file)
    scr2 <- SIGNATURE_CR
    scr2[is.na(scr2$super_target),"super_target"] <- "-"
    scr2 <- scr2[!is.element(scr2$super_target,"-"),]
    rn2 <- paste0(scr2$sample_id,"_",scr2$signature)
    rownames(scr2) <- rn2

    rn <- rn1[is.element(rn1,rn2)]
    rn <- rn[is.element(rn,rn1)]
    scr1$proper_name <- scr1$name
    scr2$proper_name <- scr2$name

    scr1$e95 <- exp(scr1$er)*qt(.975,4) - exp(scr1$er)*qt(.025,4)
    scr2$e95 <- exp(scr2$er)*qt(.975,4) - exp(scr2$er)*qt(.025,4)

    SCR1 <<- scr1[rn,]
    SCR2 <<- scr2[rn,]
  }
#    signatureConcRespPlot(SCR1[1,])
#browser()
  if(do.prep) {
    cat("  do.prep 1\n")
    scr1 <- SCR1
    scr2 <- SCR2
    res <- cbind(scr1[c("dtxsid","casrn","name","signature",
                        "sample_id","top_over_cutoff","hitcall","bmd","top","e95")],
                 scr2[c("sample_id","top_over_cutoff","hitcall","bmd","top","e95")])
    name.list <- c("dtxsid","casrn","name","signature",
                   "sid1","tc1","hc1","bmd1","top1","e951",
                   "sid2","tc2","hc2","bmd2","top2","e952")
    names(res) <- name.list

    cat("  do.prep 2\n")
    res$delta_bmd <- abs(log10(res$bmd1)-log10(res$bmd2))
    delta_tsign <- res$top1/abs(res$top1) - res$top2/abs(res$top2)
    delta_tsign <- -abs(delta_tsign)/2
    delta_tsign[delta_tsign==0] <- 1
    res$delta_top <- delta_tsign * abs(res$top1-res$top2)

    cat("  do.prep 3\n")
    res[is.na(res$tc1),"tc1"] <- 0
    res[is.na(res$tc2),"tc2"] <- 0
    res[is.na(res$hc1),"hc1"] <- 0
    res[is.na(res$hc2),"hc2"] <- 0

    cat("  do.prep 4\n")
    res$max_tc <- res$tc1
    res[res$tc2>res$tc1,"max_tc"] <- res[res$tc2>res$tc1,"tc2"]

    res$max_hc <- res$hc1
    res[res$hc2>res$hc1,"max_hc"] <- res[res$hc2>res$hc1,"hc2"]

    #res <- res[res$max_tc>1,]
    #res <- res[res$max_hc>0.5,]
    cat("  do.prep 5\n")

    RES <<- res
    cat("  do.prep 6\n")

    file <- paste0("../output/signature_replicability/method.dup.replicability.stats ",dataset,"_",sigset,".RData")
    save(RES,file=file)
    cat("  do.prep 7\n")
  }
  res <- RES
do.boxplot <- F
if(do.boxplot) {
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
}
  ####################################################################
  cat("  do cr plots 1\n")
  res <- RES

  res <- res[order(res$max_tc,decreasing=T),]
  mask <- vector(mode="integer",length=nrow(res))
  mask[] <- 0
  mask1 <- mask;mask2 <- mask
  mask1[res$tc1>1] <- 1
  mask2[res$tc2>1] <- 1
  mask12 <- mask1+mask2
  #res <- res[mask12>0,]
  cat("  do cr plots 2a:",nrow(res),"\n")

  mask <- vector(mode="integer",length=nrow(res))
  mask1 <- mask;mask2 <- mask
  mask1[res$hc1>0.5] <- 1
  mask2[res$hc2>0.5] <- 1
  mask12 <- mask1+mask2
  #res <- res[mask12>0,]
  cat("  do cr plots 2b:",nrow(res),"\n")

  mask <- vector(mode="integer",length=nrow(res))
  mask1 <- mask;mask2 <- mask
  mask1[res$bmd1<10] <- 1
  mask2[res$bmd2<10] <- 1
  mask12 <- mask1+mask2
  #res <- res[mask12>0,]
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
row.list <- unique(row.list)
  cat("  do cr plots 5:",length(row.list),"\n")

##########################################################################################
scr1 <- SCR1[row.list,]
scr2 <- SCR2[row.list,]

name.list <- c("metric","min1","min2","count")
row <- as.data.frame(matrix(nrow=1,ncol=length(name.list)))
names(row) <- name.list
mat <- NULL

##########################################################################################
hc.list <- seq(from=0,to=0.9,by=0.1)
for(hc1 in hc.list) {
	x1 <- scr1[scr1$hitcall>hc1,]
	x1 <- x1[x1$hitcall<=hc1+0.1,]
	for(hc2 in hc.list) {
	x2 <- scr2[scr2$hitcall>hc2,]
	x2 <- x2[x2$hitcall<=hc2+0.1,]
	  rn1 <- rownames(x1)
	  rn2 <- rownames(x2)
count <- length(rn1[is.element(rn1,rn2)])
row[1,"metric"] <- "hc"
row[1,"min1"] <- hc1
row[1,"min2"] <- hc2
row[1,"count"] <- count
mat <- rbind(mat,row)
}
}
##########################################################################################
tc.list <- seq(from=0,to=15,by=1)
for(tc1 in tc.list) {
	x1 <- scr1[scr1$top_over_cutoff>tc1,]
	x1 <- x1[x1$top_over_cutoff<=tc1+1,]
	for(tc2 in tc.list) {
	x2 <- scr2[scr2$top_over_cutoff>tc2,]
	x2 <- x2[x2$top_over_cutoff<=tc2+1,]
	  rn1 <- rownames(x1)
	  rn2 <- rownames(x2)
count <- length(rn1[is.element(rn1,rn2)])
row[1,"metric"] <- "tc"
row[1,"min1"] <- tc1
row[1,"min2"] <- tc2
row[1,"count"] <- count
mat <- rbind(mat,row)
}
}
##########################################################################################
ec.list <- seq(from=0,to=1.5,by=0.1)
for(ec1 in ec.list) {
	x1 <- scr1[scr1$e95>ec1,]
	x1 <- x1[x1$e95<=ec1+0.1,]
	for(ec2 in ec.list) {
	x2 <- scr2[scr2$e95>ec2,]
	x2 <- x2[x2$e95<=ec2+0.1,]
	  rn1 <- rownames(x1)
	  rn2 <- rownames(x2)
count <- length(rn1[is.element(rn1,rn2)])
row[1,"metric"] <- "e95"
row[1,"min1"] <- ec1
row[1,"min2"] <- ec2
row[1,"count"] <- count
mat <- rbind(mat,row)
}
}
##########################################################################################
bc.list <- c(1e-4,1e-3,1e-2,1e-1,1,10,100,1000)
for(bc1 in bc.list) {
	x1 <- scr1[scr1$bmd>bc1,]
	x1 <- x1[x1$bmd<=bc1*10,]
	for(bc2 in bc.list) {
	x2 <- scr2[scr2$bmd>bc2,]
	x2 <- x2[scr2$bmd<=bc2*10,]
	  rn1 <- rownames(x1)
	  rn2 <- rownames(x2)
count <- length(rn1[is.element(rn1,rn2)])
row[1,"metric"] <- "bmd"
row[1,"min1"] <- bc1
row[1,"min2"] <- bc2
row[1,"count"] <- count
mat <- rbind(mat,row)
}
}


    file <- paste0("../output/signature_replicability/method.dup.replicability.stats ",dataset,"_",sigset,".xlsx")
write.xlsx(mat,file)
if(!to.file) browser()

  par(mfrow=c(3,2),mar=c(4,4,2,2))

  plot(scr2$hitcall~scr1$hitcall,xlim=c(0,1),ylim=c(0,1),cex.lab=1.5,cex.axis=1.5,main="Hitcall",xlab="mygsea",ylab="fc",pch=".")
  lines(c(0,1),c(0,1))

  plot(scr2$top_over_cutoff~scr1$top_over_cutoff,xlim=c(0,15),ylim=c(0,15),cex.lab=1.5,cex.axis=1.5,main="Top/Cutoff",xlab="mygsea",ylab="fc",pch=".")
  lines(c(0,15),c(0,15))

  plot(scr2$e95~scr1$e95,xlim=c(0,1.5),ylim=c(0,1.5),cex.lab=1.5,cex.axis=1.5,main="Error 95% CI",xlab="mygsea",ylab="fc",pch=".")
  lines(c(0,10),c(0,10))

  plot(scr2$bmd~scr1$bmd,xlim=c(1e-3,1e3),ylim=c(1e-3,1e3),cex.lab=1.5,cex.axis=1.5,main="BMD",xlab="mygsea",ylab="fc",pch=".",log="xy")
  lines(c(0,10),c(0,10))

  if(!to.file) browser()


  par(mfrow=c(3,2),mar=c(4,4,2,2))
  cat("  cycle through signatures (rows) and run signatureConcRespPlot\n")
row.list <- row.list[1:1000]
counter <- 0
for(rn in row.list) {
    signatureConcRespPlot(scr1[rn,])
    signatureConcRespPlot(scr2[rn,])
counter <- counter+1
if(counter%%100==0) cat("finished:",counter," out of ",length(row.list),"\n")
    if(!to.file) browser()
  }  

  if(!to.file) browser()
  else dev.off()
}
