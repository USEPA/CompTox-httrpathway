#--------------------------------------------------------------------------------------
#' Look for trends in the raw data vs. the low concentraiton signature hits
#'
#' @param dir The directory from which to read all of the raw filesatalog file
#--------------------------------------------------------------------------------------
noiseHunter.rawCountsStats.vs.LowCOncSignatureHits <- function(to.file=F) {
  printCurrentFunction()

  file <- "../output/noiseHunter/noiseHunter.signatureHits.vs.rawCountsStats.xlsx"
  mat <- read.xlsx(file)
  if(to.file) {
    fname <- paste0("../output/noiseHunter/noiseHunter raw stats vs low conc hits.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,3,3,2))

  col.list <- names(mat)[11:ncol(mat)]
  for(col in col.list) {
    cat(col,"\n")
    y <- mat$nhit.0.1
    x <- mat[,col]
    pg <- mat$pg_id
    res <- lm(y ~ x + pg)
    p.col <- summary(res)[[4]][2,4]
    p.pg <- summary(res)[[4]][3,4]
    main <- paste(col,"p(raw metric,pg)",format(p.col,digits=2),format(p.pg,digits=2))
    plot(y~x,main=main,xlab=col,ylab="Sig hits <0.1 uM",cex.lab=1.2,cex.axis=1.2,ylim=c(0,2000),pch=".")
    temp <- mat[mat$pg==1,]
    y <- temp$nhit.0.1
    x <- temp[,col]
    points(x,y,pch=21,bg="cyan")
    temp <- mat[mat$pg==7,]
    y <- temp$nhit.0.1
    x <- temp[,col]
    points(x,y,pch=21,bg="red")
    temp <- mat[mat$pg==42,]
    y <- temp$nhit.0.1
    x <- temp[,col]
    points(x,y,pch=21,bg="blue")
    if(!to.file) browser()
  }

  do.by.pg <- F
  par(mfrow=c(4,2),mar=c(4,3,3,2))
  if(do.by.pg) {
    npg <- 48
    col.list <- names(mat)[11:ncol(mat)]
    for(col in col.list) {
      cat(col,"\n")
      for(pg in 1:npg) {
        temp <- mat[mat$pg_id==pg,]
        y <- temp$nhit.0.1
        x <- temp[,col]
        main <- paste("pg:",pg,col)
        plot(y~x,main=main,xlab=col,ylab="Sig hits <0.1 uM",cex.lab=1.2,cex.axis=1.2,ylim=c(0,2000))
        if(!to.file) browser()
      }
    }
  }
  if(to.file) dev.off()
}
