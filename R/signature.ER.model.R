#--------------------------------------------------------------------------------------
#' Plot the ER model results vs. eh ER pathway BMDs
#'
#' @param to.file If TRUE, send plots to a pdf
#' @param do.load If TRUE, load the data into a global
#' @param dataset The name of the data set to analyze
#' @param sigset The signature set to use
#' @param sigcatalog The signature catalog to use
#' @param method The signature scoring method
#' @export
#--------------------------------------------------------------------------------------
signature.ER.model <- function(to.file=F,
                               do.load=F,
                               sigset="pilot_large_all_100CMAP",
                               sigcatalog="signatureDB_master_catalog 2020-03-12",
                               dataset="DMEM_6hr_pilot_normal_pe_1",
                               method="mygsea"){
  printCurrentFunction()

  if(do.load) {
    annotations <<- signatureCatalogLoader(sigset,sigcatalog)
    file <- paste("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData",sep="")
    print(file)
    load(file=file)
    signaturecr <<- SIGNATURE_CR
    file <- "../input/ER/ER_model_data httr pilot.xlsx"
    ER <<- read.xlsx(file)
  }

  file <- "../input/chemicals/HTTR Pilot chemical annotation.xlsx"
  cmap <- read.xlsx(file)
  cmap <- cmap[cmap$target_key=="estrogen",]
  cmap <- cmap[is.element(cmap$dtxsid,ER$dtxsid),]
  if(to.file) {
    fname <- paste0("../input/ER/signature.ER.model.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  sig.set <- annotations[is.element(annotations$super_target,"estrogen"),"parent"]

  sigcr <- signaturecr[is.element(signaturecr$dtxsid,cmap$dtxsid),]
  sigcr <- sigcr[is.element(sigcr$signature,sig.set),c("dtxsid","signature","bmd")]

  sigcr2 <- reshape2::dcast(sigcr,signature~dtxsid,value.var="bmd")
  hitmat <- sigcr2
  hitmat[] <- 1
  hitmat[is.na(sigcr2)] <- 0
  rs <- rowSums(hitmat)
  mask <- rs
  mask[] <- 1
  mask[rs<ncol(hitmat)] <- 0
  sigcr3 <- sigcr2[mask==1,]

  res <- sigcr3
  rownames(res) <- res[,1]
  res <- res[,2:ncol(res)]
  res <- as.matrix(res)
  res <- log10(res)
  rownames(cmap) <- cmap$dtxsid
  rn <- cmap[colnames(res),"name"]
  result <- heatmap.2(as.matrix(res),
                      margins=c(15,20),
                      scale="none",
                      main="ER Signatures",
                      xlab="",
                      ylab="",
                      cexCol=1,
                      cexRow=0.3,
                      Rowv=T,
                      Colv=T,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(9,"PRGn"),
                      key.title="Key",
                      key.xlab="log10(bmd)",
                      cex.main=1,
                      labCol=rn)

  if(!to.file) browser()
  par(mfrow=c(4,2),mar=c(4,4,2,2))

  plot(c(1,1),type="n",xlim=c(1e-1,1e1),ylim=c(1e-3,1e2),log="xy",xlab="ER Model AC50",ylab="HTTR BMD",cex.lab=1.2,cex.axis=1.2)
  lines(c(1e-10,1e10),c(1e-10,1e10))
  lines(c(1e-11,1e9),c(1e-10,1e10))
  lines(c(1e-9,1e11),c(1e-10,1e10))
  for(i in 1:nrow(cmap)) {
    name <- cmap[i,"name"]
    tag <- ""
    if(name=="4-Cumylphenol") tag <- "4C"
    if(name=="4-Hydroxytamoxifen") tag <- "4H"
    if(name=="4-Nonylphenol, branched") tag <- "4NP"
    if(name=="Bisphenol A") tag <- "BPA"
    if(name=="Bisphenol B") tag <- "BPB"
    if(name=="Clomiphene citrate (1:1)") tag <- "C"
    if(name=="Fulvestrant") tag <- "F"

    dtxsid <- cmap[i,"dtxsid"]
    ac50 <- 10**(ER[is.element(ER$dtxsid,dtxsid),"AC50meanFromAUCcal"])
    temp <- sigcr3[,dtxsid]
    mval <- median(log10(temp))
    madval <- mad(log10(temp))
    bmd <- 10**mval
    bmdu <- 10**(mval+madval)
    bmdl <- 10**(mval-madval)
    cat(name,ac50,bmdu,bmd,bmdl,tag,"\n")
    lines(c(ac50,ac50),c(bmdu,bmdl))
    points(ac50,bmd,pc=21,bg="black")
    text(ac50,bmd,tag,pos=4)
  }

  if(!to.file) browser()
  else dev.off()
}
