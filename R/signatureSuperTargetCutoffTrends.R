#--------------------------------------------------------------------------------------
#' analyze the trands in PODs as a function of cutoffs
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#' mcf7_ph1_pe1_normal_all_pg
#' u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#'
#--------------------------------------------------------------------------------------
signatureSuperTargetCutoffTrends <- function(to.file=F,
                                             do.load=F,
                                             dataset="heparg2d_toxcast_pfas_pe1_normal",
                                             sigset="screen_large",
                                             method="fc",
                                             celltype="HepaRG",
                                             tcset=c(1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0),
                                             hcset=c(0.6,0.8,0.9)) {
  printCurrentFunction(paste(dataset,sigset,method))

  if(do.load) {
    res = NULL
    for(tccut in tcset) {
      for(hccut in hcset) {
        file <- paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,".xlsx")
        temp = read.xlsx(file)
        temp$hccut = hccut
        temp$tccut = tccut
        res = rbind(res,temp)
      }
    }
    RES <<- res
  }

  if(to.file) {
    fname <- paste0("../output/super_target_boxplot/",celltype,"/signatureSuperTargetCutoffTrends",celltype,"_",dataset,"_",sigset,"_",method,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,1),mar=c(4,4,2,2))
  res = RES

  chem.list = sort(unique(res$name))
  for(chem in chem.list) {
    temp1 = res[is.element(res$name,chem),]
    for(sid in unique(temp1$sid)) {
      temp2 = temp1[is.element(temp1$sid,sid),]
      plot(c(1,1),xlim=c(1,3),ylim=c(1e-3,100),log="y",xlab="T/C",ylab="POD",
           type="n",col="red",main=paste(chem,"\n",sid),cex.main=1.2)
      for(y in c(0.01,0.1,1,10,100)) lines(c(0,4),c(y,y))
      col.list=c("black","cyan","red")
      for(i in 1:length(hcset)) {
        hccut = hcset[i]
        temp3 = temp2[temp2$hccut==hccut,]
        x = temp3$tccut
        pod_st = temp3$pod_st
        pod_sig = temp3$pod_sig
        pod_bst = temp3$burst_pod_st
        pod_bsig = temp3$burst_pod_sig
        col = col.list[i]
        lines(pod_st~x,col=col,lty=1,lwd=2)
        lines(pod_sig~x,col=col,lty=2,lwd=2)

        lines(pod_bst~x,col=col,lty=1)
        lines(pod_bsig~x,col=col,lty=2)
      }
      if(!to.file) browser()
    }
  }
  if(to.file) dev.off()
}

