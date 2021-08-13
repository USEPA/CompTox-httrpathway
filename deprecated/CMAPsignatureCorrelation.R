#--------------------------------------------------------------------------------------
#' Build the correlation matrix for the CMAP signaturess
#'
#' @param min.ngene Signatures will only be saved if the number of genes is >= this value
#' @param max.ngene Signatures will only be saved if the number of genes is <= this value
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
CMAPsignatureCorrelation <- function(do.plot=F,to.file=F,ngene=100,cutoff=0.1){
  printCurrentFunction()
  load("../input/signatures/CMAP_signatures.RData")
  mat <- CMAP_signatures
  mat <- mat[mat$ngene==ngene,]
  mat.up <- mat[mat$direction=="up",]
  mat.dn <- mat[mat$direction=="dn",]
  mat.up <- mat.up[order(mat.up$signature),]
  mat.dn <- mat.dn[order(mat.dn$signature),]
  check <- mat.up$parent==mat.dn$parent
  if(sum(check) != length(check)) {
    cat("misalligenment\n")
    browser()
  }
  glist <- paste0(mat.up$gene.list,"|",mat.dn$gene.list)
  genelists = strsplit(glist, "\\|")
  print(length(genelists))
  names(genelists) = mat.up$parent

  #genelists <- genelists[1:500]
  res <- sapply(genelists,genelistCorrAll,genelists,simplify="array")
  res2 <- data.matrix(as.data.frame(res))
  res3 <- reshape2::melt(res2,id.vars=names(re2)[1],measure.vars=names(res2)[2:ncol(res2)])
  res3 <- res3[res3$value>cutoff,]
  same <- res3[,1]==res3[,2]
  res4 <- res3[!same,]
  if(do.plot) {
    if(to.file) {
      fname <- paste0("../output/signature_corr/CMAPsignatureCorrelation",ngene,".pdf")
      pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    }
    result <- heatmap.2(as.matrix(res2),
                        margins=c(10,10),
                        scale="none",
                        main="CMAP Signature Correlation ",
                        xlab="",
                        ylab="",
                        cexCol=0.7,
                        cexRow=0.7,
                        Rowv=T,
                        Colv=T,
                        trace="none",
                        hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                        key=T,
                        col=brewer.pal(9,"Reds"),
                        key.title="Key",
                        key.xlab="l2fc",
                        cex.main=1)

    if(!to.file) browser()
    else dev.off()
  }
  file <- paste0("../output/signature_corr/CMAPsignatureCorrelation",ngene,".xlsx")
  write.xlsx(res4,file)
  browser()
}
genelistCorrAll <- function(g1,genelists) {
  res <- lapply(genelists,genelistCorr,g1)
  return(res)
}

genelistCorr <- function(g1,g2) {
  g12 <- unique(c(g1,g2))
  n12 <- length(g1[is.element(g1,g2)])
  corr <- n12 / length(g12)
  return(corr)
}
