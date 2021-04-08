#--------------------------------------------------------------------------------------
#'
#' Build a heatmap of the genes in signatures in a super_target
#' @param to.file If TRUE, write plots to a file
#--------------------------------------------------------------------------------------
superTargetGeneHM <- function(to.file=F,
                              super_target="T Cell",
                              sigset="screen_large",
                              sigcatalog="signatureDB_master_catalog 2021-03-05") {
  printCurrentFunction()

  file = paste0("../input/signatures/",sigcatalog,".xlsx")

  catalog = read.xlsx(file)
  catalog = catalog[catalog[,sigset]==1,]
  catalog = catalog[is.element(catalog$super_target,super_target),]
  signature.list <- sort(catalog$signature)

  if(to.file) {
    fname <- paste0("../output/super_target_boxplot/superTargetGeneHM_",super_target,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }

  file <- paste0("../input/signatures/signatureDB_genelists.RData")
  print(file)
  load(file=file)
  # genelists

  x = NULL
  y = NULL
  #signature.list= signature.list[1:10]
  for(signature in signature.list) {
    genes = genelists[signature][[1]]
    xx = genes
    xx[] = signature
    x= c(x,xx)
    y = c(y,genes)
  }

  mat = as.data.frame(as.matrix(cbind(x,y)),stringsAsFactors=F)
  names(mat) = c("super_target","gene")
  mat$value = 1
  res= reshape2::dcast(mat,gene~super_target,fill=0)
  rownames(res) = res[,1]
  res = res[,2:ncol(res)]

  result <- heatmap.2(as.matrix(res),
                      margins=c(5,5),
                      scale="none",
                      main=paste(super_target),
                      xlab="",
                      ylab="",
                      cexCol=0.3,
                      cexRow=0.3,
                      Rowv=T,
                      Colv=T,
                      dendrogram="both",
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      key.title="Key",
                      key.xlab="-",
                      cex.main=1)

  if(!to.file) browser()

  if(to.file) dev.off()
}

