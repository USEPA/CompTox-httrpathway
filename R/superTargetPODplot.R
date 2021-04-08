#--------------------------------------------------------------------------------------
#' Generate chemicalwise boxplot of the BMD distributions by super_target
#'
#'  heparg2d_toxcast_pfas_pe1_normal
#'  mcf7_ph1_pe1_normal_block_123
#'  u2os_toxcast_pfas_pe1_normal
#'  PFAS_HepaRG
#'  PFAS_U2OS
#'  u2os_pilot_pe1_normal_null_pilot_lowconc
#'  u2os_toxcast_pfas_pe1_normal_refchems
#'  heparg2d_toxcast_pfas_pe1_normal_refchems
#'
#--------------------------------------------------------------------------------------
superTargetPODplot <- function(to.file=F,
                               dataset="PFAS_U2OS",
                               sigset="screen_large",
                               method="fc",
                               celltype="U2OS",
                               hccut=0.95,
                               tccut=1.5) {
  printCurrentFunction(paste(dataset,sigset,method))
  file = paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,"_summary.xlsx")
  mat = read.xlsx(file)
  if(to.file) {
    fname <- paste0("../output/super_target_boxplot/",celltype,"/super_target_podplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,".pdf")
    pdf(file=fname,width=5,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(2,1),mar=c(5,4,4,2))
  x = mat$burst_pod_st
  y = mat$pod_st
  plot(y~x,log="xy",xlab="POD median (uM)",ylab="POD min (uM)",xlim=c(0.0001,100),ylim=c(0.0001,100),
       cex.lab=1.2,cex.axis=1.2,main=celltype,pch=21,bg="black",cex=0.4,type="n")
  n1 = 0
  n0 = 0
  for(i in 1:nrow(mat)) {
    cex=0.4
    xx = x[i]
    yy = y[i]
    if(!is.na(xx) && !is.na(yy)) {
      col = "gray"
      if(yy/xx < 0.1) {
        col="black"
        n1 = n1+1
        stlist = mat[i,"specific_targets"]
        stlist = str_split(stlist,"\\|")[[1]]
        if(is.element("Estrogen",stlist)) {col="red";cex=0.7}
        else if(is.element("RAR",stlist)) {col="cyan";cex=0.7}
        else if(is.element("Glucocorticoid",stlist)) {col="yellow";cex=0.7}
        else if(is.element("Androgen",stlist)) {col="green";cex=0.7}
        else print(stlist)
      }
      else n0 = n0+1
      points(xx,yy,pch=21,bg=col,cex=cex)
    }
  }
  text(1e-4,50,paste("Non-specific:",n0," : Specific",n1),pos=4)
  lines(c(1e-5,1e5),c(1e-5,1e5))
  lines(c(1e-4,1e6),c(1e-5,1e5),lty=2)
  lines(c(1e-6,1e4),c(1e-5,1e5),lty=2)

  if(to.file) dev.off()
  else browser()

}

