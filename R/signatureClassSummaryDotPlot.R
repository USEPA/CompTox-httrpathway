#--------------------------------------------------------------------------------------
#'
#' Build summary plots by signature class
#' @param to.file If TRUE, write plots to a file
#' @param do.load If TRUE, load hte HTTr data file
#' @param dataset Name of the HTTr data set
#' @param sigset Name of the signature set
#' @param method Name of the scoring method
#' @param sigcatalog  Name of the signature catalog
#' @param hith.threshold Exclude rows from the dataset with hitcall less than this value
#--------------------------------------------------------------------------------------
signatureClassSummaryDotPlot <- function(to.file=F,
                                         do.load=F,
                                         dataset="DMEM_6hr_screen_normal_pe_1",
                                         sigset="screen_large",
                                         method = "mygsea",
                                         sigcatalog="signatureDB_master_catalog 2020-0505",
                                         hitcall.threshold=0.5) {
  printCurrentFunction()

  if(do.load) {
    if(dataset=="DMEM_6hr_screen_normal_pe_1") {
      file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits hits.RData")
      print(file)
      load(file)
      SIGNATURE_CR <<- SIGNATURE_CR_ACTIVE
    }
    else {
      file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
      print(file)
      load(file=file)
      SIGNATURE_CR <<- SIGNATURE_CR
    }
  }
  mat <- SIGNATURE_CR
  chems <- unique(mat[,c("dtxsid","casrn","name")])
  chems <- chems[order(chems$name),]
  mat <- mat[mat$hitcall>=hitcall.threshold,]
  mat <- mat[!is.na(mat$bmd),]

  if(to.file) {
    fname <- paste0("../output/signature_class_summary_plots/signatureClassSummaryDotPlot ",dataset,"_",sigset,"_",method,"_0.05_conthits.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,2,2))

  for(i in 1:nrow(chems)) {
    dtxsid <- chems[i,"dtxsid"]
    name <- chems[i,"name"]
    name <- chems[is.element(chems$dtxsid,dtxsid),"name"]

    temp <- mat[is.element(mat$dtxsid,dtxsid),]
    x <- temp$bmd
    top <- temp$top
    tc <- temp$top_over_cutoff
    y <- top
    y <- (top/abs(top))*tc
    #y <- tc
    plot(y~x,xlim=c(0.001,100),ylim=c(-5,5),
         xlab="BMD",ylab="Top/cutoff",log="x",
         cex.lab=1.2,cex.axis=1.2,
         main=name,type="n")
    rect(1e-10,-1,1e10,1,col="lightgray")
    lines(c(1e-10,1e10),c(0,0))
    lines(c(1e-10,1e10),c(-1,-1))
    lines(c(1e-10,1e10),c(1,1))
    lines(c(1e-10,1e10),c(-2,-2))
    lines(c(1e-10,1e10),c(2,2))
    points(y~x,pch=21,cex=0.3,bg="black")
    res <- density(y,adjust=0.2)
    yval <- res$x
    xval <- res$y + 1e-3
    lines(yval~xval,col="red")
    if(!to.file) browser()
  }
  if(to.file) dev.off()
}

