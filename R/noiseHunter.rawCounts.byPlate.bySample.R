#--------------------------------------------------------------------------------------
#' Look for trends in the raw data by plate within plate group
#'
#' @param dir The directory from which to read all of the raw filesatalog file
#--------------------------------------------------------------------------------------
noiseHunter.rawCounts.byPlate.bySample <- function(to.file=F,
                                          do.read=F,
                                          dir="../input/rawdata/mcf7_screen/") {
  printCurrentFunction()
  if(to.file) {
    fname <- paste0("../output/noiseHunter/noiseHunter.rawCounts.bySample.byPlate.mcf7.bySample.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(4,2),mar=c(4,3,3,2))

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
  smap <- smap[is.element(smap$stype,"test sample"),]

  for(pg in 1:npg) {
    cat(pg,"\n")
    temp <- res[[pg]]
    smap.pg <- smap[smap$pg_id==pg,]
    temp <- as.matrix(t(temp))
    slist <- rownames(temp)
    plate.list <- unique(smap.pg[,"plate_id"])
    for(plate in plate.list) {
      smap.pg.plate <- smap.pg[is.element(smap.pg$plate_id,plate),]
      smap.pg.plate <- smap.pg.plate[,c("sample_id","plate_id","well_id","trt_name","stype",
                                        "chem_id","dose_level","block_id","pg_id","culture_id",
                                        "n_reads","n_reads_mapd","mapd_frac","bad_probe_count",
                                        "n_cov5","n_sig80","top10_prop","gini_coef")]

      ptemp <- temp[smap.pg.plate$sample_id,]
      ptemp[ptemp>0] <- NA
      ptemp[!is.na(ptemp)] <- 1
      zc <- rowSums(ptemp,na.rm=T)
      smap.pg.plate$zero_count <- zc

      col.list <- c("n_reads","n_reads_mapd","mapd_frac","bad_probe_count",
                    "n_cov5","n_sig80","top10_prop","gini_coef")
      for(col in col.list) {
        x <- smap.pg.plate$zero_count
        y <- smap.pg.plate[,col]
        ymax <- 1.2*max(y)
        if(col=="n_sig80") ymax <- 2500
        if(col=="n_cov5") ymax <- 15000
        plot(y~x,main=paste(pg,plate,col),xlab="Zero count",ylab=col,
             xlim=c(6000,14000),ylim=c(0,ymax),cex.lab=1.2,cex.axis=1.2,pch=".")

        dx <- density(x)
        xval <- dx$x
        yscale <- 0.3*max(y)/max(dx$y)
        yval <- dx$y * yscale
        lines(yval~xval,col="red")

        dy <- density(y)
        yval <- dy$x
        yscale <- 0.3*max(x)/max(dy$y)
        xval <- 6000 + dy$y * yscale
        lines(yval~xval,col="red")
      }
      if(!to.file) browser()
    }
  }

  if(to.file) dev.off()
}
