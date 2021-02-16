#--------------------------------------------------------------------------------------
#' Compare the signature activity for a pair of chemicals
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#'
#--------------------------------------------------------------------------------------
compareTwoChems <- function(to.file=F,
                            do.load=F,
                            dataset="PFAS_U2OS",
                            sigset="screen_large",
                            method="fc",
                            celltype="U2OS",
                            dtxsid.1="DTXSID3047558",
                            dtxsid.2="DTXSID60379901",
                            hccut=0.5,
                            tccut=1) {
  printCurrentFunction(paste(dataset,sigset,method))
  if(to.file) {
    fname <- paste0("../output/signature_compare_chems/",celltype," ",dtxsid.1," ",dtxsid.2," ",dataset," ",sigset," ",method,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(2,1),mar=c(4,4,2,2))

  if(do.load) {
    file = "../output/signature_conc_resp_summary/SIGNATURE_CR_screen_large_u2os_toxcast_pfas_pe1_normal_fc_0.05_conthits.RData"
    print(file)
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    #file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_", pval,nametag,".RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    MAT <<- mat
  }
  mat = MAT
  mat.1 = mat[is.element(mat$dtxsid,dtxsid.1),]
  mat.2 = mat[is.element(mat$dtxsid,dtxsid.2),]
  mat.1 = mat.1[order(mat.1$signature),]
  mat.2 = mat.2[order(mat.2$signature),]

  mat.1[mat.1$top_over_cutoff<tccut,"bmd"] = 1000
  mat.2[mat.2$top_over_cutoff<tccut,"bmd"] = 1000

  mat.1[mat.1$hitcall<hccut,"bmd"] = 1000
  mat.2[mat.2$hitcall<hccut,"bmd"] = 1000
  x = mat.1$bmd
  x[x>1000] = 1000
  x[x<0.001] = 0.001
  y = mat.2$bmd
  y[y>1000] = 1000
  y[y<0.001] = 0.001
  name.1 = mat.1[1,"name"]
  name.2 = mat.2[1,"name"]
  plot(y~x,log="xy",xlab=dtxsid.1,ylab=dtxsid.2,main="Signature BMD",
       xlim=c(0.001,1000),ylim=c(0.001,1000),pch=21,bg="black",cex=0.5,cex.lab=1.2,cex.axis=1.2)
  lines(c(1e-4,1e4),c(1e-4,1e4))

  if(to.file) dev.off()
  else browser()

}

