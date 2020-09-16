#--------------------------------------------------------------------------------------
#'
#' Compare PODs from different conditions
#'
#' @param bmd.mode percent or abs
#' @param to.file If TRUE, write plots to a file
#--------------------------------------------------------------------------------------
signaturePodComp <- function(to.file=F,
                             sigset="screen_large",
                             method="fc") {
  printCurrentFunction()
  #' heparg2d_toxcast_pfas_pe1_normal
  #' mcf7_ph1_pe1_normal_good_pg
  #'  u2os_toxcast_pfas_pe1_normal

  if(to.file) {
    fname <- paste0("../output/super_target_boxplot/signaturePodComp_",sigset,"_",method,".pdf")
    pdf(file=fname,width=7,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,4,4))

  dataset = "heparg2d_toxcast_pfas_pe1_normal"
  celltype ="HepaRG"
  file <- paste0("../output/super_target_boxplot/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,".xlsx")
  print(file)
  heparg = read.xlsx(file)
  mat = heparg
  plot(mat$pod~mat$burst_pod,xlab="Burst POD",ylab="POD",xlim=c(0.001,100),ylim=c(0.001,100),log="xy",cex.lab=1.2,cex.axis=1.2,main=celltype,cex=0.5,pch=21,bg="white")
  lines(c(0.00001, 100000), c(0.00001,100000))
  lines(c(0.000001,10000),  c(0.00001,100000))
  lines(c(0.0001,  1000000),c(0.00001,100000))

  dataset = "u2os_toxcast_pfas_pe1_normal"
  celltype ="U2OS"
  file <- paste0("../output/super_target_boxplot/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,".xlsx")
  u2os = read.xlsx(file)
  mat = u2os
  plot(mat$pod~mat$burst_pod,xlab="Burst POD",ylab="POD",xlim=c(0.001,100),ylim=c(0.001,100),log="xy",cex.lab=1.2,cex.axis=1.2,main=celltype,cex=0.5,pch=21,bg="white")
  lines(c(0.00001, 100000), c(0.00001,100000))
  lines(c(0.000001,10000),  c(0.00001,100000))
  lines(c(0.0001,  1000000),c(0.00001,100000))

  dataset = "mcf7_ph1_pe1_normal_good_pg"
  celltype ="MCF7"
  file <- paste0("../output/super_target_boxplot/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,".xlsx")
  mcf7 = read.xlsx(file)
  mat = mcf7
  plot(mat$pod~mat$burst_pod,xlab="Burst POD",ylab="POD",xlim=c(0.001,100),ylim=c(0.001,100),log="xy",cex.lab=1.2,cex.axis=1.2,main=celltype,cex=0.5,pch=21,bg="white")
  lines(c(0.00001, 100000), c(0.00001,100000))
  lines(c(0.000001,10000),  c(0.00001,100000))
  lines(c(0.0001,  1000000),c(0.00001,100000))


  mat1 = heparg
  mat2 = u2os
  dtxsid.list = mat1$dtxsid
  dtxsid.list = dtxsid.list[is.element(dtxsid.list,mat2$dtxsid)]
  plot(c(1,1),type="n",xlab="HepaRG",ylab="U2OS",xlim=c(0.001,100),ylim=c(0.001,100),log="xy",cex.lab=1.2,cex.axis=1.2,main="POD Comparison")
  for(dtxsid in dtxsid.list) {
    x = min(mat1[is.element(mat1$dtxsid,dtxsid),"pod"])
    y = min(mat2[is.element(mat2$dtxsid,dtxsid),"pod"])
    points(x,y,cex=0.5,pch=21,bg="white")
  }
  lines(c(0.00001, 100000), c(0.00001,100000))
  lines(c(0.000001,10000),  c(0.00001,100000))
  lines(c(0.0001,  1000000),c(0.00001,100000))

  if(!to.file) browser()

  par(mfrow=c(3,2),mar=c(4,4,4,4))

  file <- "../toxcast/toxcast_pod.xlsx"
  toxcast = read.xlsx(file)

  mat1 = toxcast
  mat2 = heparg
  dtxsid.list = mat1$dtxsid
  dtxsid.list = dtxsid.list[is.element(dtxsid.list,mat2$dtxsid)]
  plot(c(1,1),type="n",xlab="ToxCast",ylab="HepaRG",xlim=c(0.001,100),ylim=c(0.001,100),log="xy",cex.lab=1.2,cex.axis=1.2,main="POD Comparison")
  for(dtxsid in dtxsid.list) {
    x = min(mat1[is.element(mat1$dtxsid,dtxsid),"pod_uM"])
    y = min(mat2[is.element(mat2$dtxsid,dtxsid),"pod"])
    points(x,y,cex=0.5,pch=21,bg="white")
  }
  lines(c(0.00001, 100000), c(0.00001,100000))
  lines(c(0.000001,10000),  c(0.00001,100000))
  lines(c(0.0001,  1000000),c(0.00001,100000))

  mat1 = toxcast
  mat2 = mcf7
  dtxsid.list = mat1$dtxsid
  dtxsid.list = dtxsid.list[is.element(dtxsid.list,mat2$dtxsid)]
  plot(c(1,1),type="n",xlab="ToxCast",ylab="MCF7",xlim=c(0.001,100),ylim=c(0.001,100),log="xy",cex.lab=1.2,cex.axis=1.2,main="POD Comparison")
  for(dtxsid in dtxsid.list) {
    x = min(mat1[is.element(mat1$dtxsid,dtxsid),"pod_uM"])
    y = min(mat2[is.element(mat2$dtxsid,dtxsid),"pod"])
    points(x,y,cex=0.5,pch=21,bg="white")
  }
  lines(c(0.00001, 100000), c(0.00001,100000))
  lines(c(0.000001,10000),  c(0.00001,100000))
  lines(c(0.0001,  1000000),c(0.00001,100000))


  mat1 = toxcast
  mat2 = u2os
  dtxsid.list = mat1$dtxsid
  dtxsid.list = dtxsid.list[is.element(dtxsid.list,mat2$dtxsid)]
  plot(c(1,1),type="n",xlab="ToxCast",ylab="U2OS",xlim=c(0.001,100),ylim=c(0.001,100),log="xy",cex.lab=1.2,cex.axis=1.2,main="POD Comparison")
  for(dtxsid in dtxsid.list) {
    x = min(mat1[is.element(mat1$dtxsid,dtxsid),"pod_uM"])
    y = min(mat2[is.element(mat2$dtxsid,dtxsid),"pod"])
    points(x,y,cex=0.5,pch=21,bg="white")
  }
  lines(c(0.00001, 100000), c(0.00001,100000))
  lines(c(0.000001,10000),  c(0.00001,100000))
  lines(c(0.0001,  1000000),c(0.00001,100000))

  if(!to.file) browser()
  if(to.file) dev.off()
}

