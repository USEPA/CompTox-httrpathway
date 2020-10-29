#--------------------------------------------------------------------------------------
#' Make plots of the data from noiseHunter1
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
noiseHunter1.plots <- function(to.file=F,
                               do.load=F,
                               dataset="DMEM_6hr_screen_normal_pe_1",
                               sigset="pilot_large_all_CMAP",
                               sigcatalog="signatureDB_master_catalog 2020-03-10",
                               method="mygsea",
                               cutoff=0.25){
  printCurrentFunction()

  if(do.load) {
    file <- paste("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData",sep="")
    print(file)
    load(file=file)
    signaturecr <<- SIGNATURE_CR
    annotations <<- signatureCatalogLoader(sigset,sigcatalog)

    file <- paste0("../input/chemicals/",dataset,"_chemical_map.xlsx")
    cmap <<- read.xlsx(file)
  }
  if(to.file) {
    fname <- paste0("../output/noiseHunter/noiseHunter1 ",cutoff,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(4,1),mar=c(4,4,2,2))

  file <- paste0("../output/noiseHunter/noiseHunter1.chemicals ",cutoff,".xlsx")
  res <- read.xlsx(file)
  x <- res$pg_id
  y <- res$nsig
  boxplot(y~x,xlab="Plate group",main="Number of low conc. hits per chemical by plate group",cex.lab=1.2,cex.axis=0.5,ylim=c(0,500))

  x <- unique(res$pg_id)
  y <- x
  y[] <- 0
  for(i in 1:length(x)) {
    pg_id <- x[i]
    temp <- res[res$pg_id==pg_id,]
    y[i] <- nrow(temp)
  }
  barplot(y~x,xlab="Plate group",main="Number of chemicals with low conc hits by plate group",cex.lab=1.2,cex.axis=1.2,cex.names=0.5)

  x <- res$block_id
  y <- res$nsig
  boxplot(y~x,xlab="Block",main="Number of low conc. hits per chemical by block",cex.lab=1.2,cex.axis=1.2,ylim=c(0,500))

  x <- unique(res$block_id)
  y <- x
  y[] <- 0
  for(i in 1:length(x)) {
    block_id <- x[i]
    temp <- res[res$block_id==block_id,]
    y[i] <- nrow(temp)
  }
  barplot(y~x,xlab="Block",main="Number of chemicals with low conc hits by block",cex.lab=1.2,cex.axis=1.2,cex.names=1.2)

  file <- paste0("../output/noiseHunter/noiseHunter1.signatures ",cutoff,".xlsx")
  res <- read.xlsx(file)

  if(!to.file) browser()

  if(to.file) dev.off()
}
