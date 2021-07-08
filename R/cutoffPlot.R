#####################################################################################################
#' Calculate the signature-wise cutoffs based on the analytical method
#' which does not break any correlations between genes
#'
#' @param basedir Directory that holds FCMAT2 and CHEM_DICT files.
#' @param dataset Name of actual dataset to base cutoff on.
#' @param sigcatalog The name of the signature catalog to use
#' @param sigset THe signature set
#' @param method The scoring method, either fc or gsea
#' @param pval The p-value for the baseline distribution
#' @param seed Random seed.
#' @param nlowconc Only include the lowest nlowconc concentrations for each chemical
#' @param mc.cores NUmber of coresto use when running parallel
#' @param dtxsid.exclude dtxsids to exclude, default NULL
#' @param do.load If TRUE, reload the FCMAT2 matrix, signature catalog and chemical dictionary, and store in globals
#' @param do.cov If TRUE, calculate the covariance matrix and store in a global
#' @param do.compare If TRUE, compare the cutoffs with those from the original method with no gene-gene correlation
#' @param to.file If TRUE, and do.compare=TRUE, send a plot of the comparison to a file
#' @param verbose If TRUE, write a line for each signature to show progress.
#'
#' @return No output.
#' @export
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
