#--------------------------------------------------------------------------------------
#' Plot summary data for the PFAS immune  data
#'
#'  heparg2d_toxcast_pfas_pe1_normal
#'  u2os_toxcast_pfas_pe1_normal
#'
#'  PFAS_HepaRG
#'  PFAS_U2OS
#'
#--------------------------------------------------------------------------------------
pfasImmuneSummaryPlot <- function(to.file=F,
                                  do.load=F,
                                  sigcatalog="signatureDB_master_catalog 2021-03-05",
                                  dataset.other="heparg2d_toxcast_pfas_pe1_normal",
                                  dataset.pfas="PFAS_HepaRG",
                                  sigset="screen_large",
                                  method="fc",
                                  celltype="HepaRG",
                                  bmdcut=1000,
                                  hccut=0.9,
                                  tccut=1.5) {
  printCurrentFunction(paste(celltype,sigset))


  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset.other,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat.other = SIGNATURE_CR
     file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset.pfas,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat.pfas = SIGNATURE_CR
    dlist = unique(mat.pfas$dtxsid)
    mat.other = mat.other[!is.element(mat.other$dtxsid,dlist),]

    mat.other = mat.other[mat.other$top_over_cutoff>1,]
    mat.pfas = mat.pfas[mat.pfas$top_over_cutoff>1,]

    MAT.PFAS <<- mat.pfas
    MAT.OTHER <<- mat.other
  }
  mat.other = MAT.OTHER
  mat.pfas = MAT.PFAS
  mat.other$cclass = "Other"
  mat.pfas$cclass = "PFAS"
  file = paste0("../input/signatures/",sigcatalog,".xlsx")
  catalog = read.xlsx(file)
  if(to.file) {
    fname <- paste0("../output/PFAS/pfasImmuneSummary_",celltype,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(5,4,4,2))
  st.list = sort(unique(catalog[is.element(catalog$target_class,"Immune"),"super_target"]))

  name.list = c("sample_id","dtxsid","name","super_target","signature","bmd","top_over_cutoff","top","hitcall","cclass")
  row = as.data.frame(matrix(nrow=1,ncol=length(name.list)))
  names(row) = name.list
  res = NULL
  for(st in st.list) {
    temp.other = mat.other[is.element(mat.other$super_target,st),]
    temp.pfas = mat.pfas[is.element(mat.pfas$super_target,st),]
    cat(st,nrow(temp.pfas),nrow(temp.other),"\n")

    x = temp.other[temp.other$top_over_cutoff>tccut,name.list]
    #cat(nrow(x),"\n")
    x = x[x$bmd<bmdcut,]
    #cat(nrow(x),"\n")
    x = x[x$hitcall>hccut,]
    #cat(nrow(x),"\n")
    if(nrow(x)>0) res = rbind(res,x)
    x = temp.pfas[temp.pfas$top_over_cutoff>tccut,name.list]
    #cat(nrow(x),"\n")
    x = x[x$bmd<bmdcut,]
    #cat(nrow(x),"\n")
    x = x[x$hitcall>hccut,]
    #cat(nrow(x),"\n")
    if(nrow(x)>0) res = rbind(res,x)
    #if(!is.null(res)) browser()
    if(nrow(temp.other)>0 && nrow(temp.pfas)>0) {

      x = log10(temp.other$bmd)
      y = temp.other$top_over_cutoff
      y[y>20] = 20

      qx = quantile(x,probs=seq(0,1,0.05),na.rm=T)
      qy = quantile(y,probs=seq(0,1,0.05),na.rm=T)
      xcut = qx[2]
      ycut = qy[20]

      xcut = 1
      ycut = 2.5
      xmask = x
      xmask[] = 0
      xmask[x<xcut] = 1
      ymask = y
      ymask[] = 0
      ymask[y>ycut] = 1
      n1 = sum(xmask*ymask) / length(x)
      n2 = sum(xmask * (1-ymask)) / length(x)
      n3 = sum((1-xmask) * ymask) / length(x)
      n4 = sum((1-xmask) * (1-ymask)) / length(x)

      plot(y~x,xlab="log(bmd uM)",ylab="T/C",main=paste(st,"\nnon-PFAS"),cex.lab=1.2,cex.axis=1.2,xlim=c(-4,2),ylim=c(1,20),pch=21,cex=0.2)
      lines(c(xcut,xcut),c(-10,100))
      lines(c(-10,10),c(ycut,ycut))
      text(-4,15,format(n1,digits=2),pos=4)
      text(-4,1,format(n2,digits=2),pos=4)
      text(xcut,15,format(n3,digits=2),pos=4)

      x = log10(temp.pfas$bmd)
      y = temp.pfas$top_over_cutoff
      y[y>20] = 20
      plot(y~x,xlab="log(bmd uM)",ylab="T/C",main=paste(st,"\nPFAS"),cex.lab=1.2,cex.axis=1.2,xlim=c(-4,2),ylim=c(1,20),pch=21,cex=0.2)
      lines(c(xcut,xcut),c(-10,100))
      lines(c(-10,10),c(ycut,ycut))
      xmask = x
      xmask[] = 0
      xmask[x<xcut] = 1
      ymask = y
      ymask[] = 0
      ymask[y>ycut] = 1
      n1 = sum(xmask*ymask) / length(x)
      n2 = sum(xmask * (1-ymask)) / length(x)
      n3 = sum((1-xmask) * ymask) / length(x)
      n4 = sum((1-xmask) * (1-ymask)) / length(x)
      text(-4,15,format(n1,digits=2),pos=4)
      text(-4,1,format(n2,digits=2),pos=4)
      text(xcut,15,format(n3,digits=2),pos=4)
      if(!to.file) browser()
    }
  }
  if(to.file) dev.off()
  file = paste0("../output/PFAS/pfasImmuneSummary_specific_",celltype,".xlsx")
  write.xlsx(res,file)

}
