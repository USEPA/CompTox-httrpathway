#####################################################################################################
#' Calculate the signature-wise cutoffs based on the analytical method
#' which does not break any correlations between genes
#'
#' @param dataset Name of actual dataset to base cutoff on.
#' @param sigset THe signature set
#' @param method The scoring method, either fc or gsea
#' @param pval The p-value for the baseline distribution
#' @param nlowconc Only include the lowest nlowconc concentrations for each chemical
#' @param to.file If TRUE, and do.compare=TRUE, send a plot of the comparison to a file
#' @importFrom graphics par plot lines hist
#' @importFrom grDevices pdf dev.off
#' @importFrom openxlsx read.xlsx
#'
#' @return No output.
#' @export cutoffPlot
#####################################################################################################
cutoffPlot = function(to.file=F,
                      dataset="mcf7_ph1_pe1_normal_block_123_allPG",
                      sigset="screen_large",
                      method="fc",
                      pval=0.05,
                      nlowconc=2){
  printCurrentFunction()
  if(to.file) {
    file = paste0("../output/signature_cutoff/signature_cutoff_",sigset,"_",dataset,"_",method,"_",pval,"_",nlowconc,"_with_gene_correlations.pdf")
    pdf(file=file,width=10,height=8,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(2,2),mar=c(5,4,4,3))

  file = paste0("../output/signature_cutoff/signature_cutoff_",sigset,"_",dataset,"_",method,"_",pval,".xlsx")
  print(file)
  oldcuts = read.xlsx(file)

  file = paste0("../output/signature_cutoff/signature_cutoff_",sigset,"_",dataset,"_",method,"_",pval,"_",nlowconc,"_with_gene_correlations.xlsx")
  print(file)
  newcuts = read.xlsx(file)

  rownames(oldcuts) = oldcuts$signature
  rownames(newcuts) = newcuts$signature
  siglist = rownames(oldcuts)
  siglist = siglist[is.element(siglist,rownames(newcuts))]
  oldcuts = oldcuts[siglist,]
  newcuts = newcuts[siglist,]
  x = oldcuts$cutoff
  y = newcuts$cutoff
  plot(y~x,xlim=c(0,0.25),ylim=c(0,0.25),xlab="old method",ylab="new method",pch=".",cex.lab=1.2,cex.lab=1.2,main=paste(dataset,nlowconc))
  lines(c(0,0.25),c(0,0.25))
  delta = y-x
  breaks = seq(from=-0.4,to=0.4,by=0.02)
  hist(delta,breaks=breaks,cex.lab=1.2,cex.lab=1.2,xlab="Delta",main=paste(dataset,nlowconc))
  if(!to.file) browser()
  else dev.off()
}
