#--------------------------------------------------------------------------------------
#' Make plots of the data by plate
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
noiseHunter1.pgplots.raw.v3 <- function(to.file=F){
  printCurrentFunction()
  if(to.file) {
    fname <- paste0("../output/noiseHunter/noiseHunter raw pg distribution by plate.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(5,4),mar=c(4,3,3,2))
  basedir="../input/httr_mcf7_screen/rawfc/platewise/"
  for(pg_id in c(4,7)) {
    cat("Plate group",pg_id,"\n")
    dataset <- paste0("mcf7_screen_raw_l2fc_pg_",pg_id)

    file <- paste0(basedir,"pg",pg_id,"_plateLevel_FCMAT2.RData")
    print(file)
    load(file)
    data.name <- paste0("pg",pg_id,"_plate_dat_wide2")
    FCMAT2 <- get(data.name)
    fcmat2 <- FCMAT2
    fcmat2[is.na(fcmat2)] <- 0
    rn <- rownames(fcmat2)
    temp <- str_split(rn,"_")
    temp2 <- as.data.frame(do.call(rbind,temp),stringsAsFactors=F)
    temp2$rn <- rn
    pid.list <- unique(temp2[,1])
    for(pid in pid.list) {
      sk <- temp2[is.element(temp2[,1],pid),"rn"]
      fcmat_pid <- fcmat2[sk,]
      temp <- fcmat_pid
      cmean <- colMeans(temp)
      cmed <- apply(temp,FUN=median,MARGIN=2)
      sk <- skewness(cmean)
      plot(density(cmean),xlim=c(-2,2),xlab="mean(gene-wise l2fc)",ylim=c(0,3),cex.lab=0.8,cex.axis=0.8,cex.main=0.8,
           main=paste("PG mean",pg_id,pid,"\nSK=",format(sk,digits=2),"SD=",format(sd(cmean),digits=2)))
      plot(density(cmed),xlim=c(-1,1),xlab="median(gene-wise l2fc)",ylim=c(0,3),cex.lab=0.8,cex.axis=0.8,cex.main=0.8,
           main=paste("PG median",pg_id,pid))
    }
  }
  if(!to.file) browser()

  if(to.file) dev.off()
}
