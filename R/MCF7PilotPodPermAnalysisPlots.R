#' Carry out analyses of different POD methods for the pilot study
#' Do hte permutation analyses
#'
#'
#' * MCF7_pilot_DMEM_6hr_pilot_normal_pe_1
#' * MCF7_pilot_DMEM_12hr_pilot_normal_pe_1
#' * MCF7_pilot_DMEM_24hr_pilot_normal_pe_1
#' * MCF7_pilot_PRF_6hr_pilot_normal_pe_1
#' * MCF7_pilot_PRF_12hr_pilot_normal_pe_1
#' * MCF7_pilot_PRF_24hr_pilot_normal_pe_1

MCF7PilotPodPermAnalysisPlots <- function(to.file=F,
                                          dataset="MCF7_pilot_DMEM_6hr_pilot_normal_pe_1",
                                          method="gsea",
                                          celltype="MCF7",
                                          sigset="screen_large",
                                          hccut=0.9,
                                          tccut=1,
                                          cutoff=3,
                                          pval=0.05) {
  printCurrentFunction()
  dir = "../output/mcf7_pilot/"
  file = paste0(dir,"/signature_perm_pod_",sigset,"_",dataset,"_",method,"_",hccut,"_",cutoff,".xlsx")
  mat = read.xlsx(file)
  mat = mat[order(mat$name),]
  dlist = unique(mat$dtxsid)

  if(to.file) {
    fname = paste0(dir,"MCF7PilotPodPermAnalysisPlots_",sigset,"_",dataset,"_",method,"_",hccut,"_",cutoff,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(4,4),mar=c(4,3,3,2))

  for(dtxsid in dlist) {
    temp = mat[is.element(mat$dtxsid,dtxsid),]
    x = temp$nsigtot
    y = temp$signature_pod_95
    plot(y~x,xlab="Signatures",ylab="POD",log="y",xlim=c(0,15000),ylim=c(0.001,1000),main=paste(temp[1,"name"],"\n",dataset),pch=21,cex=0.2,cex.main=0.8)
    if(!to.file) browser()
  }
  if(to.file) dev.off()
}
