#--------------------------------------------------------------------------------------
#' Look at quadrants of test plates to see if there is a checkerboard pattern
#' in the read counts
#'
#' @param to.file If TRUE, write plots to a file
#-------------------------------------------------------------------------------------
httr.plate.quad.reads <- function(to.file=F,set="mcf7") {
  printCurrentFunction()
  if(to.file) {
    fname <- paste0("../httr/",set," httr.plate.quad.reads.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(4,3),mar=c(4,4,2,2))

  if(set=="mcf7") {
    file <- "../httr/httr_mcf7_ph1_metadata.xlsx"
    ymin <- 19
    ymax <- 23
    npg <- 48
  }
  if(set=="u2os") {
    file <- "../httr/httr_u2os_pfas _metadata.xlsx"
    ymin <- 19
    ymax <- 23
    npg <- 4
  }
  mat <- read.xlsx(file)
  row.odd <-  c("A","C","E","G","I","K","M","O")
  row.even <- c("B","D","F","H","J","L","N","P")
  col.even <- seq(from=1,to=24,by=2)
  col.odd <- seq(from=2,to=24,by=2)
  q1 <- NULL
  for(row in row.odd) {
    for(col in col.even) {
      scol <- as.character(col)
      if(nchar(scol)==1) scol <- paste0("0",scol)
      q1 <- c(q1,paste0(row,scol))
    }
  }
  q2 <- NULL
  for(row in row.odd) {
    for(col in col.odd) {
      scol <- as.character(col)
      if(nchar(scol)==1) scol <- paste0("0",scol)
      q2 <- c(q2,paste0(row,scol))
    }
  }
  q3 <- NULL
  for(row in row.even) {
    for(col in col.odd) {
      scol <- as.character(col)
      if(nchar(scol)==1) scol <- paste0("0",scol)
      q3 <- c(q3,paste0(row,scol))
    }
  }
  q4 <- NULL
  for(row in row.even) {
    for(col in col.even) {
      scol <- as.character(col)
      if(nchar(scol)==1) scol <- paste0("0",scol)
      q4 <- c(q4,paste0(row,scol))
    }
  }


  for(pg in 1:npg) {
    temp <- mat[mat$pg_id==pg,]
    block <- temp[1,"block_id"]
    plate.list <- unique(temp$plate_id)
    for(plate in plate.list) {
      ptemp <- temp[is.element(temp$plate_id,plate),]

      x <- NULL
      y <- NULL

      dmso <- ptemp[is.element(ptemp$stype,"vehicle control"),]
      yy <- dmso$n_reads_mapd
      xx <- yy
      xx[] <- "dmso"
      x <- c(x,xx)
      y <- c(y,yy)

      ptemp <- ptemp[is.element(ptemp$stype,"test sample"),]
      ptemp <- ptemp[order(ptemp$well_id),]
      yy <- ptemp[is.element(ptemp$well_id,q1),"n_reads_mapd"]
      xx <- yy
      xx[] <- "Q1"
      x <- c(x,xx)
      y <- c(y,yy)
      yy <- ptemp[is.element(ptemp$well_id,q2),"n_reads_mapd"]
      xx <- yy
      xx[] <- "Q2"
      x <- c(x,xx)
      y <- c(y,yy)
      yy <- ptemp[is.element(ptemp$well_id,q3),"n_reads_mapd"]
      xx <- yy
      xx[] <- "Q3"
      x <- c(x,xx)
      y <- c(y,yy)
      yy <- ptemp[is.element(ptemp$well_id,q4),"n_reads_mapd"]
      xx <- yy
      xx[] <- "Q4"
      x <- c(x,xx)
      y <- c(y,yy)

      y <- log2(y)
      boxplot(y~x,main=paste("B,PG,P",block,pg,plate),ylim=c(ymin,ymax),
              cex.axis=1.2,cex.lab=1.2,ylab="n_reads_mapd",xlab="Quadrant")
    }
    if(!to.file) browser()
  }
  if(to.file) dev.off()
}
