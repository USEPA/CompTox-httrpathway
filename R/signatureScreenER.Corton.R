#--------------------------------------------------------------------------------------
#'
#' Plot the ER model vs. the Corton ER signatures
#--------------------------------------------------------------------------------------
signatureScreenER.Corton <- function(to.file=F,
                                     do.load=F,
                                     do.prep=F,
                                     dataset="mcf7_ph1_pe1_normal_block_123",
                                     sigset="stress_test",
                                     sigcatalog="signatureDB_master_catalog 2021-03-02",
                                     method="fc") {
  printCurrentFunction()
  if(to.file) {
    fname <- paste0("../output/signature_ER/signature_ER_",dataset,"_",sigset,"_",method,"_Corton.pdf")
    pdf(file=fname,width=5,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(2,1),mar=c(5,4,4,3))

  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    mat = mat[mat$hitcall>0.9,]
    mat = mat[mat$top_over_cutoff>1.5,]

    SIGNATURE_CR <<- mat
  }

  if(do.prep) {
    screen <- SIGNATURE_CR
    file <- "../output/signature_ER/ER_model_data 2020-03-26.xlsx"
    file <- "../input/ER/ER_model_data 2019-11-01.xlsx"
    erdata <- read.xlsx(file)
    erdata <- erdata[erdata$AC50meanFromAUC<2,]
    dtxsid.list <- erdata$dtxsid
    dtxsid.list <- dtxsid.list[is.element(dtxsid.list,screen$dtxsid)]
    erdata <- erdata[is.element(erdata$dtxsid,dtxsid.list),]
    name.list <- c("dtxsid","casrn","name","AC50meanFromAUC")
    mat <- erdata[,name.list]
    mat$ermodel.AC50 <- mat$AC50meanFromAUC

    screen <- screen[is.element(screen$super_target,"Estrogen"),]
    screen <- screen[is.element(screen$dtxsid,dtxsid.list),]

    x <- screen[is.element(screen$dtxsid,dtxsid.list[1:10]),]
    x <- x[x$bmd<1,]
    y <- as.data.frame(table(x$signature))
    y <- y[order(y$Freq,decreasing=T),]
    sig.list <- y[y$Freq>=3,"Var1"]
    screen <- screen[is.element(screen$signature,sig.list),]
    mat$ERsignature.median <- NA
    mat$ERsignature.uci <- NA
    mat$ERsignature.lci <- NA

    for(i in 1:nrow(mat)) {
      dtxsid <- mat[i,"dtxsid"]
      bmd <- median(screen[is.element(screen$dtxsid,dtxsid),"bmd"])
      lbmd <- log10(bmd)
      cat(dtxsid,lbmd,"\n")
      if(!is.na(lbmd)) {
        q <- quantile(lbmd,probs=seq(0,1,0.025))
        mat[i,"ERsignature.lci"] <- q[2]
        mat[i,"ERsignature.uci"] <- q[40]
        mat[i,"ERsignature.median"] <- q[21]
      }
    }
    mat = mat[!is.na(mat$ERsignature.median),]
    file <- paste0("../output/signature_ER/signature_ER_",dataset,"_",sigset,"_",method,".xlsx")
    write.xlsx(mat,file)
    MAT.ER <<- mat
  }

  mat <- MAT.ER

  x <- mat[,"AC50meanFromAUC"]
  y <- mat[,"ERsignature.median"]
  res <- lm(y~x)
  sres <- summary(res)
  r2 <- sres$adj.r.squared
  rmse <- sres$sigma
  main <- paste0("R2=",format(r2,digits=2),"   RMSE=",format(rmse,digits=2))

  plot(c(1,1),type="n",xlab="ER Model log(AC50, uM)",ylab="HTTr log(BMD, uM)",
       xlim=c(-3,2),ylim=c(-3,2),cex.axis=1.2,cex.lab=1.2,
       main=main)
  lines(c(-5,5),c(-5,5))
  for(i in 1:nrow(mat)) {
    x <- mat[i,"AC50meanFromAUC"]
    y <- mat[i,"ERsignature.median"]
    uci <- mat[i,"ERsignature.uci"]
    lci <- mat[i,"ERsignature.lci"]
    lines(c(x,x),c(uci,lci))
    points(x,y,pch=21,cex=0.6,bg="black")
  }

  if(to.file) dev.off()
  else browser()
}
