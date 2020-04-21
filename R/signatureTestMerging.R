#--------------------------------------------------------------------------------------
#'
#' Examine the merging code effect
#' @param to.file If TRUE, write plots to a file
#--------------------------------------------------------------------------------------
signatureTestMerging <- function(to.file=F,
                        dataset="DMEM_6hr_pilot_normal_pe_1",
                        sigset="pilot_tiny",
                        method="mygsea",
                        cname="Fulvestrant") {
  printCurrentFunction()

  if(to.file) {
    fname <- paste0("../output/signature_score_summary/signatureTestMerging",dataset,"_",sigset," ",method," ",cname,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,2,2))

  file <- paste0("../output/signature_score_summary/signaturescoremat_",sigset,"_",dataset,"_",method,".RData")
  print(file)
  load(file=file)
  smat1 <- signaturescoremat

  file <- paste0("../output/signature_score_summary/signaturescoremat_",sigset,"_",dataset,"_",method,"_directional.RData")
  print(file)
  load(file=file)
  smat0 <- signaturescoremat

  smat1 <- smat1[is.element(smat1$name,cname),]
  smat0 <- smat0[is.element(smat0$name,cname),]

  par(mfrow=c(3,2),mar=c(4,4,2,2))
  for(signature in smat1$signature) {
    if(substr(signature,1,4)=="CMAP") {
      smat1both <- smat1[is.element(smat1$signature,signature),]
      sig <- paste0(signature," up")
      smat0up <- smat0[is.element(smat0$signature,sig),]
      sig <- paste0(signature," dn")
      smat0dn <- smat0[is.element(smat0$signature,sig),]

      x <- smat1both$conc
      y <- smat1both$signature_score
      plot(y~x,xlab="Conc",log="x",xlim=c(0.01,100),ylim=c(-1,1),main=paste(cname,"\n",signature),
           cex.lab=1.2,cex.axis=1.2,ylab="MyGSEA Score",pch=21,bg="black",cex=1.5)
      lines(c(1e-5,1e5),c(0,0))
      x <- smat0up$conc#*1.1
      y <- smat0up$signature_score
      points(y~x,pch=24,bg="green")
      x <- smat0dn$conc#*0.9
      y <- smat0dn$signature_score
      points(y~x,pch=25,bg="red")
      if(!to.file) browser()
    }
  }
  if(to.file) dev.off()
}

