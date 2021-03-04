#' Replicate Chemical signature Plot
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#' u2os_pilot_pe1_normal_null_pilot_lowconc
#' u2os_toxcast_pfas_pe1_normal_refchems
#' heparg2d_toxcast_pfas_pe1_normal_refchems
signatureReplication <- function(to.file=F,
                                 dataset="mcf7_ph1_pe1_normal_block_123",
                                 sigset="screen_large",
                                 method="fc",
                                 celltype="MCF7",
                                 hccut=0.95,
                                 tccut=1.5) {
  printCurrentFunction(paste(dataset,sigset,method))

  file = paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,"_summary.xlsx")
  mat = read.xlsx(file)
  if(to.file) {
    fname <- paste0("../output/super_target_boxplot/",celltype,"/signature_replication_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,".pdf")
    pdf(file=fname,width=5,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(2,1),mar=c(5,4,4,2))

  x = log10(mat$pod_st)
  x = x[x<3]
  plot(density(x),xlim=c(-3,2),cex.axis=1.3,cex.lab=1.4,xlab="log(bmd uM)",main="All Chemicals")
  dlist = mat$dtxsid
  dups = dlist[duplicated(dlist)]
  name.list = c("dtxsid","s1","s2")
  res = as.data.frame(matrix(nrow=length(dups),ncol=length(name.list)))
  names(res) = name.list
  res$dtxsid=dups
  for(i in 1:nrow(res)) {
    dtxsid = res[i,"dtxsid"]
    temp = mat[is.element(mat$dtxsid,dtxsid),]
    res[i,"s1"] = temp[1,"pod_st"]
    res[i,"s2"] = temp[2,"pod_st"]
  }
  x = log10(res$s1)
  y = log10(res$s2)
  plot(y~x,xlim=c(-1,3),ylim=c(-1,3),cex.axis=1.3,cex.lab=1.4,
       xlab="log(bmd uM, rep 1)",ylab="log(bmd uM, rep 2)",main="Replicates",
       type = "n")
  for(i in 1:length(y)) {
    xx = x[i]
    yy = y[i]
    col = "gray"
    if(yy<3 && xx<3) col="green"
    else if(yy==3 && xx<3) col="orange"
    else if(yy<3 && xx==3) col="orange"
    else col="blue"
    points(xx,yy,pch=21,bg=col,cex=0.8)
  }
  lines(c(-5,5),c(-5,5))
  lines(c(-6,4),c(-5,5),lty=2)
  lines(c(-4,6),c(-5,5),lty=2)

  if(!to.file) browser()
  else dev.off()

}
