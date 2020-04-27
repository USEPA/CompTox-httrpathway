#--------------------------------------------------------------------------------------
#' Compare the MCF7 biomodality data between raw counts and l2fc
#'
#' @param to.file If TRUE, send plots to a pdf
#' @param do.load If TRUE, load the data into a global
#' @param dataset The name of the data set to analyze
#' @param sigset The signature set to use
#' @param sigcatalog The signature catalog to use
#' @param method The signature scoring method
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
noiseHunter.l2fc.raw.stats <- function(to.file=F,dataset="DMEM_6hr_screen_normal_pe_1"){
  printCurrentFunction()
  file <- paste0("../output/noiseHunter/noiseHunter ",dataset," pg distribution.xlsx")
  data.l2fc <- read.xlsx(file)
  file <- paste0("../output/noiseHunter/noiseHunter.rawCounts.bySample mcf7.xlsx")
  data.raw <- read.xlsx(file)
  if(to.file) {
    fname <- paste0("../output/noiseHunter/noiseHunter.l2fc.raw.stats.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(5,4,4,2))

  names.l2fc <- names(data.l2fc)[2:ncol(data.l2fc)]
  names.raw <- names(data.raw)[2:ncol(data.raw)]
  for(nraw in names.raw) {
    x <- data.raw[,nraw]
    for(nl2fc in names.l2fc) {
      y <- data.l2fc[,nl2fc]
      stats <- lm(y~x)
      xx <- summary(stats)
      pval <- xx[[4]][2,4]

      plot(y~x,xlab=nraw,ylab=nl2fc,main=paste(nl2fc,"~",nraw,":",format(pval,digits=2)),cex.lab=1.2,cex.axis=1.2,type="n")
       for(i in 1:length(x)) {
        color <- "black"
        if(is.element(i,c(1,2,7,42,45))) color <- "red"
        text(x[i],y[i],i,col=color)
      }
      if(!to.file) browser()
    }
  }
  if(to.file) dev.off()
}
