#' Export the Null set thresholds
#'
#' Null dataset and dataset should have already been scored using signatureScore
#' and the given sigset and method.
#' @param sigset Name of the signature set.
#' @param dataset Name of the data set.
#' @param method Pathway scoring method in c("fc", "gsva", "mygsea")
#' @param pval Desired cutoff p-value.
#'
mcf7_all_good_pg_null_threshold <- function(to.file=F,
                                            sigset="screen_large",
                                            dataset1="mcf7_ph1_pe1_normal_all_pg",
                                            dataset2="mcf7_ph1_pe1_normal_good_pg",
                                            method="fc") {
  printCurrentFunction()

  file <- paste0("../output/signature_score_summary/signaturescoremat_",sigset,"_",dataset1,"_",method,"_nullSetThreshold.RData")
  load(file)
  all_pg = signature_null_thresholds
  file <- paste0("../output/signature_score_summary/signaturescoremat_",sigset,"_",dataset2,"_",method,"_nullSetThreshold.RData")
  load(file)
  good_pg = signature_null_thresholds

  if(to.file) {
    fname <- paste0("../output/signature_score_summary/mcf7_all_good_pg_null_threshold.pdf")
    pdf(file=fname,width=5.5,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(2,1),mar=c(4,5,2,2))

  sigs = unique(all_pg$signature)
  sigs = sigs[is.element(sigs,good_pg$signature)]

  #################################################################
  all_pg_pval = all_pg[is.element(all_pg$type,"pval"),]
  rownames(all_pg_pval) = all_pg_pval$signature
  all_pg_pval = all_pg_pval[sigs,]

  good_pg_pval = good_pg[is.element(good_pg$type,"pval"),]
  rownames(good_pg_pval) = good_pg_pval$signature
  good_pg_pval = good_pg_pval[sigs,]

  plot(main="pval",good_pg_pval$cutoff ~ all_pg_pval$cutoff,xlab="All PG",ylab="Good PG",cex.axis=1.2,cex.lab=1.2,xlim=c(0.01,0.6),ylim=c(0.01,0.6),log=("xy"))
  meandiff = mean(good_pg_pval$cutoff/all_pg_pval$cutoff,na.rm=T)
  text(0.01,0.5,paste("Mean(good/all)=",format(meandiff,digits=2)),pos=4)
  dall = density(all_pg_pval$cutoff,na.rm=T)
  dgood = density(good_pg_pval$cutoff,na.rm=T)
  lines(c(0.001,1),c(0.001,1))
  x = dall$x
  y = 0.01*dall$y+0.01
  lines(y~x)
  x = 0.01*dgood$y + 0.01
  y =dgood$x
  lines(y~x)

  #################################################################
  all_pg_pval = all_pg[is.element(all_pg$type,"numsd"),]
  rownames(all_pg_pval) = all_pg_pval$signature
  all_pg_pval = all_pg_pval[sigs,]

  good_pg_pval = good_pg[is.element(good_pg$type,"numsd"),]
  rownames(good_pg_pval) = good_pg_pval$signature
  good_pg_pval = good_pg_pval[sigs,]

  plot(main="numsd",good_pg_pval$cutoff ~ all_pg_pval$cutoff,xlab="All PG",ylab="Good PG",cex.axis=1.2,cex.lab=1.2,xlim=c(0,0.3),ylim=c(0,0.3))
  dall = density(all_pg_pval$cutoff,na.rm=T)
  dgood = density(good_pg_pval$cutoff,na.rm=T)
  lines(c(0,1),c(0,1))
  lines(0.01*dall$y~dall$x)
  x = 0.01*dgood$y
  y =dgood$x
  lines(y~x)


    if(!to.file) browser()
  else dev.off()
}

