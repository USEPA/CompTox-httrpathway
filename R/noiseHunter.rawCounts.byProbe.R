#--------------------------------------------------------------------------------------
#' Look for trends in the raw data by plate group and probe
#'
#' @param dir The directory from which to read all of the raw filesatalog file
#--------------------------------------------------------------------------------------
noiseHunter.rawCounts.byProbe <- function(to.file=F,
                                          do.read=F,
                                          step1=F,
                                          dir="../input/rawdata/mcf7_screen/") {
  printCurrentFunction()
  if(to.file) {
    fname <- paste0("../output/noiseHunter/noiseHunter.rawCounts.byProbe mcf7.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(5,4),mar=c(3,2,2,2))

  npg <- 48
  if(do.read) {
    res <- list()
    for(pg in 1:npg) {
      cat("pg:",pg,"\n")
      file <- paste0(dir,"counts_pg_",pg,".csv")
      temp <- read.csv(file)
      rownames(temp) <- temp[,1]
      temp <- temp[,2:ncol(temp)]
      temp$pg <- pg
      res[[pg]] <- temp
    }
    RES <<- res
  }
  res <- RES

  if(step1) {
    result <- NULL
    for(pg in 1:npg) {
      cat(pg,"\n")
      temp <- res[[pg]]
      temp <- temp[,2:(ncol(temp)-1)]

      temp <- as.matrix(t(temp))
      temp[temp>0] <- NA
      temp[!is.na(temp)] <- 1
      cs <- colSums(temp,na.rm=T)
      result <- rbind(result,cs)
    }
    RESULT1 <<- result
  }
  result1 <- RESULT1
  plist <- colnames(result1)
  nprobe <- ncol(result1)
  name.list <- c("probe_id","mean.in","mean.out","delta","pvalue")
  result2 <- as.data.frame(matrix(nrow=nprobe,ncol=length(name.list)))
  names(result2) <- name.list
  in.list <- c(1,2,6,9,7,17,42,45)
  out.list <- seq(from=1,to=48,by=1)
  out.list <- out.list[!is.element(out.list,in.list)]
  for(i in 1:length(plist)) {
    pid <- plist[i]
    x <- result1[,pid]
    xin <- x[in.list]
    xout <- x[out.list]
    result2[i,"probe_id"] <- pid
    result2[i,"mean.in"] <- mean(xin)
    result2[i,"mean.out"] <- mean(xout)
    result2[i,"delta"] <- (mean(xin)-mean(xout))/mean(xout)
    stats <- wilcox.test(xin,xout,alternative="greater",exact=F)
    p <- stats$p.value
    result2[i,"pvalue"] <- p
    if(p<0.0001) {
      dx <- density(x,adjust=0.01)
      ymax <- max(dx$y)
      plot(dx,main=paste(pid,"\n",i,":",format(p,digits=2)),xlim=c(0,max(dx$x)))
      for(j in in.list) lines(c(x[j],x[j]),c(0,ymax),col="red")
    if(!to.file) browser()
    }
    if(i%%1000==0) cat("finished ",i," out of ",length(plist),"\n")
  }
  file <- paste0("../output/noiseHunter/noiseHunter.rawCounts.byProbe mcf7.xlsx")
  write.xlsx(result2,file)
  if(to.file) dev.off()
}
