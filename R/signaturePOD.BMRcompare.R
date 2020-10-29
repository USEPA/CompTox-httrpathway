#--------------------------------------------------------------------------------------
#'
#' Compare teh PODs with different BMR values
#' @param to.file If TRUE, write plots to a file
#' @param dataset Name of data set.
#' @param sigset Name of signature set.
#' @param method Pathway scoring method in c("fc", "gsva", "mygsea")
#' @param bmr_scale	bmr scaling factor. Default = 1.349

#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#' u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#--------------------------------------------------------------------------------------
signaturePOD.BMRcompare <- function(to.file=F,
                                    dataset="mcf7_ph1_pe1_normal_block_123",
                                    sigset="screen_large",
                                    method="fc",
                                    bmr_scale=1,
                                    hccut=0.9) {
  printCurrentFunction()

  if(to.file) {
    fname <- paste0("../output/signature_pod/signaturePOD.BMRcompare ",dataset,"_",sigset,"_",method,".pdf")
   pdf(file=fname,width=5,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(2,1),mar=c(4,4,2,2))
  file <- paste0("../output/signature_pod/signature_pod_",sigset,"_",dataset,"_",method,"_",hccut,".xlsx")
  mat.bmr.0 = read.xlsx(file)
  file <- paste0("../output/signature_pod/signature_pod_",sigset,"_",dataset,"_",method,"_bmr_scale_",bmr_scale,"_",hccut,".xlsx")
  mat.bmr = read.xlsx(file)

  x = mat.bmr.0$signature_pod_95
  y = mat.bmr$signature_pod_95
  ratio=x/y
  main=paste(dataset,"\nRatio=",format(mean(ratio),digits=2))
  plot(y~x,xlab="BMR=1.349",ylab=paste0("BMR=",bmr_scale),log="xy",pch=21,cex=0.5,xlim=c(0.01,100),ylim=c(0.01,100),main=main)
  lines(c(0.001,1000),c(0.001,1000))
  lines(c(0.0001,100),c(0.001,1000))
  lines(c(0.01,10000),c(0.001,1000))
  if(!to.file) browser()
  else dev.off()
}

