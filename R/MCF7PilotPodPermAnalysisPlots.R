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

  file = paste0(dir,"ER_chems all mcf7_ph1_pe1_normal_block_123_allPG screen_large 0.9 10.xlsx")
  ermodel = read.xlsx(file)
  erchems = c(
    "Fulvestrant",
    "4-Hydroxytamoxifen",
    "Clomiphene citrate (1:1)",
    "Bisphenol B",
    "Bisphenol A",
    "4-Nonylphenol, branched",
    "4-Cumylphenol"
  )
  ermodel = ermodel[is.element(ermodel$name,erchems),]
  rownames(ermodel) = ermodel$name
  ermodel = ermodel[erchems,]
  rownames(ermodel) = ermodel$dtxsid
  nchem = nrow(ermodel)
  #dlist = dlist[is.element(dlist,ermodel$dtxsid)]
  dlist = ermodel$dtxsid

  if(to.file) {
    fname = paste0(dir,"MCF7PilotPodPermAnalysisPlots_",sigset,"_",dataset,"_",method,"_",hccut,"_",cutoff,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(4,3),mar=c(4,4,3,2))

  for(dtxsid in dlist) {
    temp = mat[is.element(mat$dtxsid,dtxsid),]
    x = temp$nsigtot
    y = temp$signature_pod_95
    plot(y~x,xlab="Signatures",ylab="uM",log="y",xlim=c(0,15000),ylim=c(0.001,1000),las=1,
         main=paste(temp[1,"name"]),pch=21,cex=0.2,cex.lab=1.2,cex.axis=1.2)
    pod = min(ermodel[dtxsid,"hts.pod.agonist"],ermodel[dtxsid,"hts.pod.antagonist"])
    pod = 10**pod
    lines(c(0,1000000),c(pod,pod),lwd=2,col="red")
  }
  if(!to.file) browser()
  if(to.file) dev.off()
}
