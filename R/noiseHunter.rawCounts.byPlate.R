#--------------------------------------------------------------------------------------
#' Look for trends in the raw data by plate within plate group
#'
#' @param dir The directory from which to read all of the raw filesatalog file
#--------------------------------------------------------------------------------------
noiseHunter.rawCounts.byPlate <- function(to.file=F,
                                          do.read=F,
                                          dir="../input/rawdata/mcf7_screen/") {
  printCurrentFunction()
  if(to.file) {
    fname <- paste0("../output/noiseHunter/noiseHunter.rawCounts.bySample.byPlate.mcf7.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(5,4,4,2))

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
    file <- "../output/noiseHunter/httr_mcf7_ph1_metadata.xlsx"
    smap <- read.xlsx(file)
    SMAP <<- smap
  }
  res <- RES
  smap <- SMAP
  rownames(smap) <- smap$sample_id
  plate.list <- unique(smap$plate_id)
  plate.list <- plate.list[!is.na(plate.list)]
  nplate <- length(plate.list)
  name.list <- c("pg","plate","zero.count.mean","zero.count.sd","nonzero.mean","nonzero.sd","readcount.mean","readcount.sd")
  result <- as.data.frame(matrix(nrow=nplate,ncol=length(name.list)))
  names(result) <- name.list
  counter <- 0
  for(pg in 1:npg) {
    cat(pg,"\n")
    temp <- res[[pg]]

    temp <- as.matrix(t(temp))
    slist <- rownames(temp)
    tmap <- smap[slist,"plate_id"]
    tmap <- tmap[!is.na(tmap)]
    plist <- unique(tmap)
    for(plate in plist) {
      counter <- counter+1

      ptemp <- temp[is.element(plist,plate),]

      temp2 <- ptemp
      temp3 <- temp2
      temp3[temp3==0] <- NA
      rm3 <- rowMeans(temp3,na.rm=T)
      temp4 <- ptemp
      temp4[temp4>0] <- NA
      temp4[!is.na(temp4)] <- 1
      rs4 <- rowSums(temp4,na.rm=T)
      rs4.mean <- mean(rs4)
      rs4.sd  <- sd(rs4)

      result[counter,"pg"] <- pg
      result[counter,"zero.count.mean"] <- rs4.mean
      result[counter,"zero.count.sd"] <- rs4.sd
      result[counter,"nonzero.mean"] <- mean(rm3)
      result[counter,"nonzero.sd"] <- sd(rm3)

      main <- paste0("pg: ",pg," : ",plate," (",format(rs4.mean,digits=1),":",format(rs4.sd,digits=2),")")
      plot(rs4~rm3,xlab="Mean counts for non-zeros",ylab="Number of zeros",
           main=main,cex.lab=1.2,cex.axis=1.2,log="x",
           xlim=c(1e1,1e3),ylim=c(0,20000),pch=".")

      x <- density(rm3)
      xval <- x$x
      yval <- x$y * 1000000
      lines(yval~xval,col="red")

      x <- density(rs4)
      yval <- x$x
      xval <- 10 + x$y * 100000
      lines(yval~xval,col="red")

      lrm3 <- log10(rm3)
      mod <- lm(rs4~lrm3)
      slope <- mod$coefficients[2]
      intercept <- mod$coefficients[1]
      x1 <- 1
      x2 <- 1000
      y1 <- intercept
      y2 <- intercept + slope*(log10(x2)-log10(x1))
      lines(c(x1,x2),c(y1,y2),col="blue")

      rs <- log10(rowSums(ptemp))
      result[counter,"readcount.mean"] <- mean(rs)
      result[counter,"readcount.sd"] <- sd(rs)

      main <- paste0("pg: ",pg," : ",plate," (",format(mean(rs),digits=2),":",format(sd(rs),digits=2),")")
      plot(rs4~rs,xlab="log10(mean counts)",ylab="Number of zeros",
           main=main,cex.lab=1.2,cex.axis=1.2,
           xlim=c(3,8),ylim=c(0,20000),pch=".")

      x <- density(rs)
      xval <- x$x
      yval <- x$y * 5000
      lines(yval~xval,col="red")

      x <- density(rs4)
      yval <- x$x
      xval <- 3 + x$y * 3000
      lines(yval~xval,col="red")

      mod <- lm(rs4~rs)
      slope <- mod$coefficients[2]
      intercept <- mod$coefficients[1]
      x1 <- 0
      x2 <- 10
      y1 <- intercept
      y2 <- intercept + slope*(x2-x1)
      lines(c(x1,x2),c(y1,y2),col="blue")

      if(!to.file) browser()
    }
  }
  file <- paste0("../output/noiseHunter/noiseHunter.rawCounts.bySample.byPlate.mcf7.xlsx")
  write.xlsx(result,file)

  if(to.file) dev.off()
}
