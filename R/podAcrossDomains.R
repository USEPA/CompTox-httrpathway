#--------------------------------------------------------------------------------------
#' Analyze the PODs across domains
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#'
#--------------------------------------------------------------------------------------
podAcrossDomains <- function(to.file=F) {
  printCurrentFunction()

  if(!exists("DSSTOX")) {
    file = "../input/DSSTox/DSSTox 2020-01-02.RData"
    load(file=file)
    rownames(DSSTOX) = DSSTOX$dtxsid
    DSSTOX <<- DSSTOX
  }

  sigset="screen_large"
  method="fc"
  hccut = 0.95
  tccut = 1.5

  dataset="mcf7_ph1_pe1_normal_block_123"
  celltype="MCF7"
  file = paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,"_summary.xlsx")
  httr.mcf7 = read.xlsx(file)

  dataset="u2os_toxcast_pfas_pe1_normal"
  celltype="U2OS"
  file = paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,"_summary.xlsx")
  httr.u2os = read.xlsx(file)

  dataset="heparg2d_toxcast_pfas_pe1_normal"
  celltype="HepaRG"
  file = paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,"_summary.xlsx")
  httr.heparg = read.xlsx(file)

  file = "../toxcast/toxcast_pod.xlsx"
  toxcast = read.xlsx(file)

  file = "../htpp/htpp_bmc.xlsx"
  htpp = read.xlsx(file)
  htpp[is.na(htpp$cv_bmc),"cv_bmc"] = 1000
  htpp.mcf7 = htpp[is.element(htpp$celltype,"MCF7"),]
  htpp.u2os = htpp[is.element(htpp$celltype,"U2OS"),]
  dtxsid.list = unique(
    c(
      httr.mcf7$dtxsid,
      httr.u2os$dtxsid,
      httr.heparg$dtxsid,
      htpp$dtxsid,
      toxcast$dtxsid
    )
  )
  cat(length(dtxsid.list),"\n")
  dtxsid.list = dtxsid.list[!is.na(dtxsid.list)]
  cat(length(dtxsid.list),"\n")

  chems = DSSTOX[dtxsid.list,c("dtxsid","casrn","name")]
  rownames(chems) = chems$dtxsid

  mat = chems
  mat$toxcast.pod = NA
  mat$httr.mcf7.sig.bmd = NA
  mat$httr.mcf7.st.bmd = NA
  mat$httr.u2os.sig.bmd = NA
  mat$httr.u2os.st.bmd = NA
  mat$httr.heparg.sig.bmd = NA
  mat$httr.heparg.st.bmd = NA
  mat$htpp.mcf7.bmc = NA
  mat$htpp.u2os.bmc = NA

  for(dtxsid in mat$dtxsid) {
    if(is.element(dtxsid,toxcast$dtxsid)) {
      mat[dtxsid,"toxcast.pod"] = toxcast[is.element(toxcast$dtxsid,dtxsid),"pod_uM"]
    }
    if(is.element(dtxsid,httr.mcf7$dtxsid)) {
      temp = httr.mcf7[is.element(httr.mcf7$dtxsid,dtxsid),]
      mat[dtxsid,"httr.mcf7.sig.bmd"] = median(temp$pod_sig)
      mat[dtxsid,"httr.mcf7.st.bmd"] = median(temp$pod_st)
    }
    if(is.element(dtxsid,httr.heparg$dtxsid)) {
      temp = httr.heparg[is.element(httr.heparg$dtxsid,dtxsid),]
      mat[dtxsid,"httr.heparg.sig.bmd"] = median(temp$pod_sig)
      mat[dtxsid,"httr.heparg.st.bmd"] = median(temp$pod_st)
    }
    if(is.element(dtxsid,httr.u2os$dtxsid)) {
      temp = httr.u2os[is.element(httr.u2os$dtxsid,dtxsid),]
      mat[dtxsid,"httr.u2os.sig.bmd"] = median(temp$pod_sig)
      mat[dtxsid,"httr.u2os.st.bmd"] = median(temp$pod_st)
    }
    if(is.element(dtxsid,htpp.u2os$dtxsid)) {
      temp = htpp.u2os[is.element(htpp.u2os$dtxsid,dtxsid),]
      mat[dtxsid,"htpp.u2os.bmc"] = median(temp$cv_bmc)
    }
    if(is.element(dtxsid,htpp.mcf7$dtxsid)) {
      temp = htpp.mcf7[is.element(htpp.mcf7$dtxsid,dtxsid),]
      mat[dtxsid,"htpp.mcf7.bmc"] = median(temp$cv_bmc)
    }
  }

  if(to.file) {
    fname = paste0("../output/super_target_boxplot/podAcrossDomains ",Sys.Date(),".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,4,2))
  color.list = c("red","orange","yellow","green","blue","violet","black","gray","white","brown")
  col.list = names(mat)[5:ncol(mat)]
  temp = mat[,col.list]
  temp[is.na(temp)] = 1000
  x = rowMin(temp)
  mat$min.nam = x
  col.list = names(mat)[5:ncol(mat)]
  ic = 0
  for(col in col.list) {
    plot(c(1,1),type="n",xlim=c(0.001,1000),ylim=c(0.001,1000),xlab="ToxCast",ylab=col,main=col,log="xy",cex.lab=1.2,cex.axis=1.2)
    lines(c(1e-5,1e5),c(1e-5,1e5))
    lines(c(1e-5,1e5),c(1e-6,1e4))
    lines(c(1e-5,1e5),c(1e-4,1e6))
    ic = ic+1
    x = mat$toxcast.pod
    y = mat[,col]
    points(y~x,pch=21,bg=color.list[ic],cex=0.5)
    if(!to.file) browser()
  }
  if(!to.file) browser()
  else dev.off()
  file <- paste0("../output/super_target_boxplot/podAcrossDomains ",Sys.Date(),".xlsx")
  write.xlsx(mat,file)
}

