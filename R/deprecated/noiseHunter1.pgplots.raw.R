#--------------------------------------------------------------------------------------
#' Make plots of the data by plate group
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
noiseHunter1.pgplots.raw <- function(to.file=F, do.swap=F){
  printCurrentFunction()
  if(to.file) {
    fname <- paste0("../output/noiseHunter/noiseHunter1 raw pg distribution.pdf")
    if(do.swap) fname <- paste0("../output/noiseHunter/noiseHunter1 raw pg distribution swap.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(5,4),mar=c(4,3,3,2))
  basedir="../input/fcdata/"
  for(pg_id in c(4,7)) {
    cat("Plate group",pg_id,"\n")
    dataset <- paste0("mcf7_screen_raw_l2fc_pg_",pg_id)

    file <- paste0(basedir,"FCMAT2_",dataset,".RData")
    print(file)
    load(file)
    fcmat2 <- FCMAT2
    fcmat2[is.na(fcmat2)] <- 0

    temp <- fcmat2
    cmean <- colMeans(temp)
    cmed <- apply(temp,FUN=median,MARGIN=2)
    sk <- skewness(cmean)
    plot(density(cmean),xlim=c(-2,2),xlab="mean(gene-wise l2fc)",ylim=c(0,5),cex.lab=0.8,cex.axis=0.8,cex.main=0.8,
         main=paste("PG mean",pg_id,"\nSK=",format(sk,digits=2),"SD=",format(sd(cmean),digits=2)))
    plot(density(cmed),xlim=c(-1,1),xlab="median(gene-wise l2fc)",ylim=c(0,5),cex.lab=0.8,cex.axis=0.8,cex.main=0.8,
         main=paste("PG median",pg_id))
  }
  if(!to.file) browser()

  if(to.file) dev.off()
}
