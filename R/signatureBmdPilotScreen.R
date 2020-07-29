#--------------------------------------------------------------------------------------
#'
#' Get the signature ranks for chemicals
#--------------------------------------------------------------------------------------
signatureBmdPilotScreen <- function(to.file=F,
                                    do.prep=F,
                                    do.load.screen=F,do.load.pilot=F,
                                    dataset.screen="DMEM_6hr_screen_normal_pe_1",
                                    dataset.pilot="DMEM_6hr_pilot_normal_pe_1",
                                    sigset.screen="screen_large",
                                    sigset.pilot="pilot_large",
                                    sigcatalog="signatureDB_master_catalog 2020-04-05",
                                    method="mygsea") {
  printCurrentFunction()
  if(to.file) {
    fname <- paste0("../output/signature_pilot_screen/signature_pilot_screen_",dataset.screen,"_",sigset.screen,"_",method,".pdf")
    pdf(file=fname,width=5,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(2,1),mar=c(5,4,4,3))

  if(do.load.screen) {
    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset.screen,"_",dataset.screen,"_",method,"_0.05_conthits hits.RData")
    print(file)
    load(file)
    SIGNATURE_CR.screen <<- SIG_CONC_RESPONSE_ACTIVE
  }

  if(do.load.pilot) {
    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset.pilot,"_",dataset.pilot,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file)
    SIGNATURE_CR.pilot <<- SIGNATURE_CR
  }

  if(do.prep) {
    screen <- SIGNATURE_CR.screen
    pilot <- SIGNATURE_CR.pilot
    pilot <- pilot[pilot$hitcall>0.5,]

    dtxsid.list <- unique(screen$dtxsid)
    dtxsid.list <- dtxsid.list[is.element(dtxsid.list,pilot$dtxsid)]
    screen <- screen[is.element(screen$dtxsid,dtxsid.list),]
    pilot <- pilot[is.element(pilot$dtxsid,dtxsid.list),]

    sig.list <- unique(screen$signature)
    sig.list <- sig.list[is.element(sig.list,pilot$signature)]
    screen <- screen[is.element(screen$signature,sig.list),]
    pilot <- pilot[is.element(pilot$signature,sig.list),]

    mat <- screen[,c("dtxsid","signature","bmd")]
    names(mat)[3] <- "bmd.screen"
    mat$bmd.pilot <- NA
    for(i in 1:nrow(mat)) {
      dtxsid <- mat[i,"dtxsid"]
      signature <- mat[i,"signature"]
      temp <- pilot[is.element(pilot$dtxsid,dtxsid),]
      temp <- temp[is.element(temp$signature,signature),]
      if(nrow(temp)>0) {
        mat[i,"bmd.pilot"] <- temp[1,"bmd"]
      }
    }
    mat <- mat[!is.na(mat$bmd.pilot),]
    BMD.MAT <<- mat
  }
  mat <- BMD.MAT
  x <- mat$bmd.pilot
  y <- mat$bmd.screen

  lx <- log10(x)
  ly <- log10(y)
  res <- lm(ly~lx)
  sres <- summary(res)
  r2 <- sres$adj.r.squared
  rmse <- sres$sigma
  main <- paste0("R2=",format(r2,digits=2),"   RMSE=",format(rmse,digits=2))
  plot(y~x,log="xy",xlab="BMD(pilot)",ylab="BMD(screen)",main=main,
       cex.lab=1.2,cex.axis=1.2,xlim=c(0.001,100),ylim=c(0.001,100),
       cex=0.2)
  lines(c(1e-6,1000),c(1e-6,1000))
  if(to.file) dev.off()
  else browser()
}
