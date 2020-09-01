#--------------------------------------------------------------------------------------
#' plot the number of extreme hits vs. number of genes in the signature
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#' u2os_toxcast_pfas_pe1_normal
#--------------------------------------------------------------------------------------
signatureChemStatsPlot <- function(to.file=F,
                                   do.prep=F,
                                   sigcatalog="signatureDB_master_catalog 2020-08-14",

                                   sigset="screen_large",
                                   method="fc") {
  printCurrentFunction()
  if(do.prep) {
    file <- paste0("../input/signatures/",sigcatalog,".xlsx")
    catalog0 <- read.xlsx(file)
    catalog0 <- catalog0[catalog0[,sigset]==1,]
    plist <- unique(catalog0$parent)
    np <- length(plist)
    catalog <- as.data.frame(matrix(nrow=np,ncol=2))
    names(catalog) <- c("signature","ngene")
    catalog$signature <- plist
    rownames(catalog) <- catalog$signature
    for(sig in plist) {
      temp <- catalog0[is.element(catalog0$parent,sig),"ngene"]
      catalog[sig,"ngene"] <- sum(temp)
    }

    celltype <- "U2OS"
    dataset <- "u2os_toxcast_pfas_pe1_normal"
    file <- paste0("../output/signature_cluster/",celltype,"/signature_stats_",celltype,"_",sigset,"_",dataset,"_",method,".xlsx")
    u2os <- read.xlsx(file)
    rownames(u2os) <- u2os$signature

    celltype <- "MCF7"
    dataset <- "mcf7_ph1_pe1_normal_good_pg"
    file <- paste0("../output/signature_cluster/",celltype,"/signature_stats_",celltype,"_",sigset,"_",dataset,"_",method,".xlsx")
    mcf7 <- read.xlsx(file)
    rownames(mcf7) <- mcf7$signature

    celltype <- "HepaRG"
    dataset <- "heparg2d_toxcast_pfas_pe1_normal"
    file <- paste0("../output/signature_cluster/",celltype,"/signature_stats_",celltype,"_",sigset,"_",dataset,"_",method,".xlsx")
    heparg <- read.xlsx(file)
    rownames(heparg) <- heparg$signature

    slist <- u2os$signature
    slist <- slist[is.element(slist,mcf7$signature)]
    slist <- slist[is.element(slist,heparg$signature)]
    slist <- slist[is.element(slist,catalog$signature)]
    u2os <- u2os[slist,]
    mcf7 <- mcf7[slist,]
    heparg <- heparg[slist,]
    catalog <- catalog[slist,]
    name.list <- c("signature","ngene","u2os","mcf7","heparg")
    res <- as.data.frame(matrix(nrow=length(slist),ncol=length(name.list)))
    names(res) <- name.list
    res[,1] <- slist
    res[,2] <- catalog$ngene
    res[,3] <- u2os$fhigh
    res[,4] <- mcf7$fhigh
    res[,5] <- heparg$fhigh
    RES <<- res
  }
  res <- RES
  if(to.file) {
    fname <- paste0("../output/signature_cluster/signature_ngene_hits_allcelltype.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(2,1),mar=c(4,4,2,2))

  plot(c(0,0),type="n",xlim=c(0,1000),ylim=c(0,0.2),
       xlab="Genes in Signature",ylab="f(chemicals active)",
       cex.lab=1.5,cex.axis=1.5)
  points(res$u2os~res$ngene,pch=21,cex=0.5,bg="red")
  points(res$mcf7~res$ngene,pch=21,cex=0.5,bg="black")
  points(res$heparg~res$ngene,pch=21,cex=0.5,bg="cyan")

  d1 <- density(res$u2os)
  y <- d1$x
  x <- d1$y*10
  lines(y~x,col="red",lwd=2)

  d1 <- density(res$heparg)
  y <- d1$x
  x <- d1$y*10
  lines(y~x,col="cyan",lwd=2)

  d1 <- density(res$mcf7)
  y <- d1$x
  x <- d1$y*10
  lines(y~x,col="black",lwd=2)

  if(!to.file) browser()
  else dev.off()

}
