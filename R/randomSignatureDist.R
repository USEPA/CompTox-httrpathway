#--------------------------------------------------------------------------------------
#' Examine the distribution of hits for the random signatures
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#' mcf7_ph1_pe1_normal_all_pg
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#'
#--------------------------------------------------------------------------------------
randomSignatureDist <- function(to.file=F,
                                do.load=F,
                                dataset="mcf7_ph1_pe1_normal_good_pg",
                                sigset="screen_large",
                                method="fc",
                                celltype="U2MCF7OS") {
  printCurrentFunction(paste(dataset,sigset,method))

  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    mat = mat[is.element(mat$super_target,"Random"),]

    MAT <<- mat
  }
  mat = MAT
  if(to.file) {
    fname <- paste0("../output/signature_conc_resp_summary/randomSignatureDist",celltype,"_",dataset,"_",sigset,"_",method,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,2,2))

  v1 = mat$top_over_cutoff
  v2 = mat$bmd
  v3 = mat$hitcall
  boxplot(v1,main=paste("T/C",celltype),cex.lab=1.2,cex.axis=1.2,xlab="T/C")
  boxplot(v2,main=paste("BMD",celltype),cex.lab=1.2,cex.axis=1.2,xlab="BMD",log="y")
  boxplot(v3,main=paste("HC",celltype),cex.lab=1.2,cex.axis=1.2,xlab="HC")

  plot(v1~v3,xlab="HC",ylab="T/C",cex.lab=1.2,cex.axis=1.2,main=celltype)
  cat("T/C\n")
  print(quantile(v1,probs=seq(0,1,0.05),na.rm=T))
  print(range(v1,na.rm=T))

  cat("HC\n")
  print(quantile(v3,probs=seq(0,1,0.05),na.rm=T))
  print(range(v3,na.rm=T))

    if(!to.file) browser()
  else dev.off()
}

