#--------------------------------------------------------------------------------------
#' Heatmap of the low conc rows
#'
#' @param min.ngene Signatures will only be saved if the number of genes is >= this value
#' @param max.ngene Signatures will only be saved if the number of genes is <= this value
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
method.rep.hm <- function(to.file=F,
                          dataset="mcf7_ph1_pe1_normal_good_pg",
                          sigset="screen_large") {
  printCurrentFunction()
  if(to.file) {
    fname <- paste0("../output/signature_replicability/method.rep.hm ",dataset,"_",sigset,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
    file <- paste0("../output/signature_replicability/method.dup.replicability.stats ",dataset,"_",sigset,".xlsx")
mat <- read.xlsx(file)
mlist <- c("hc","tc","e95","bmd")
for(metric in mlist) {
temp <- mat[is.element(mat$metric,metric),]
mat2 <- reshape2::dcast(temp,min1~min2,value.var="count")
rownames(mat2) <- mat2[,1]
mat2 <- mat2[,2:ncol(mat2)]
mat2[mat2==0] <- 1
mat2 <- log10(mat2)
    result <- heatmap.2(as.matrix(mat2),
                        margins=c(4,4),
                        dendrogram="none",
                        scale="none",
                        main=metric,
                        xlab="fc",
                        ylab="mygsea",
                        cexCol=1,
                        cexRow=1,
                        Rowv=F,
                        Colv=F,
                        trace="none",
                        hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                        key=T,
                        col=brewer.pal(9,"Reds"),
                        key.title="Key",
                        key.xlab="log10(counts)",
                        cex.main=1)

if(!to.file) browser()

}
if(to.file) dev.off()
}
