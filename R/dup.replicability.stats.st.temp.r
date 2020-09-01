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

    res2 <- RES2
    res2$et1 <- res2$tc1 / res2$ec1
    res2$et2 <- res2$tc2 / res2$ec2

    RES2 <<- res2
  }
  res2 <- RES2
  plot(res2$tc2~res2$tc1,pch=".")
  plot(res2$hc2~res2$hc1,pch=".")
  plot(res2$ec2~res2$ec1,pch=".")
  plot(res2$bmd2~res2$bmd1,log="xy",pch=".")
  plot(res2$caikwt2~res2$caikwt1,pch=".")
  plot(res2$et2~res2$et1,pch=".",xlim=c(0,6),ylim=c(0,6))

  x <- vector(length=nrow(res2),mode="integer")
  y <- x

  cat("hitcall>0.5\n")
  x[] <- 0; y[] <- 0
  x[res2$hc1>0.5] <- 1
  y[res2$hc2>0.5] <- 1
  a <- sum(x*y); b <- sum(x*(1-y)); c <- sum((1-x)*y); d <- sum((1-x)*(1-y)); txt <- TxT(a,b,c,d)
  print(txt$mat)

  cat("tc>2\n")
  x[] <- 0; y[] <- 0
  x[res2$tc1>2] <- 1
  y[res2$tc2>2] <- 1
  a <- sum(x*y); b <- sum(x*(1-y)); c <- sum((1-x)*y); d <- sum((1-x)*(1-y)); txt <- TxT(a,b,c,d)
  print(txt$mat)

  cat("et>2\n")
  x[] <- 0; y[] <- 0
  x[res2$et1>2] <- 1
  y[res2$et2>2] <- 1
  a <- sum(x*y); b <- sum(x*(1-y)); c <- sum((1-x)*y); d <- sum((1-x)*(1-y)); txt <- TxT(a,b,c,d)
  print(txt$mat)

  title <- "hc>0.9 and tc>2"
  cat(title,"\n")
  x[] <- 0; y[] <- 0
  mask1 <- x
  mask2 <- x
  mask3 <- x
  mask1[res2$tc1>3] <- 1
  mask2[res2$hc1>0.9] <- 1
  mask3[res2$caikwt1<0.1] <- 1
  x <- mask1*mask3

  mask1 <- y
  mask2 <- y
  mask3 <- y
  mask1[res2$tc2>3] <- 1
  mask2[res2$hc2>0.9] <- 1
  mask3[res2$caikwt2<0.1] <- 1
  y <- mask1*mask3
  a <- sum(x*y); b <- sum(x*(1-y)); c <- sum((1-x)*y); d <- sum((1-x)*(1-y)); txt <- TxT(a,b,c,d,rowname=title)
  print(txt$mat)

    if(!to.file) browser()
  else dev.off()
}
