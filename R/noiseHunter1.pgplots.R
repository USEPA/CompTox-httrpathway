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
noiseHunter1.pgplots <- function(to.file=F,
                                 do.load=F,
                                 dataset="DMEM_6hr_screen_normal_pe_1"){
  printCurrentFunction()
  basedir="../input/fcdata/"
  if(do.load) {
    ds2 <- str_replace(dataset,"_pgnorm","")
    file <- paste0("../input/chemicals/",ds2,"_chemical_map.xlsx")
    cmap <<- read.xlsx(file)

    file <- paste0(basedir,"FCMAT2_",dataset,".RData")
    print(file)
    load(file)
    FCMAT2 <<- FCMAT2

    file <- paste0(basedir,"FCMAT2.PV.",dataset,".RData")
    print(file)
    load(file)
    FCMAT2.PV <<- FCMAT2.PV

    file <- paste0(basedir,"FCMAT2.SE.",dataset,".RData")
    print(file)
    load(file)
    FCMAT2.SE <<- FCMAT2.SE
  }

  #---------------------------------------------------------------------------------
  # se
  #---------------------------------------------------------------------------------
  if(to.file) {
    fname <- paste0("../output/noiseHunter/noiseHunter ",dataset," pg distribution.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(5,3),mar=c(4,3,3,2))

  fcmat2 <- FCMAT2
  fcmat2[is.na(fcmat2)] <- 0
  fcmat2.se <- FCMAT2.SE
  fcmat2.se[is.na(fcmat2.se)] <- 0
  fcmat2.pv <- FCMAT2.PV
  fcmat2.pv[is.na(fcmat2.pv)] <- 1

  gene.list <- colnames(fcmat2)

  rownames(cmap) <- cmap$sample_key
  pg.list <- sort(unique(cmap$pg_id))

  ngene <- ncol(FCMAT2)
  output <- NULL
  for(i in 1:length(pg.list)) {
    pg_id <- pg.list[i]
    cat("Plate group",pg_id,"\n")
    sk.list <- cmap[cmap$pg_id==pg_id,"sample_key"]

    temp <- fcmat2[sk.list,]
    cmean <- colMeans(temp)
    sk <- skewness(cmean)
    plot(density(cmean),xlim=c(-1,1),xlab="mean(gene-wise l2fc)",ylim=c(0,15),cex.lab=0.8,cex.axis=0.8,cex.main=0.8,
         main=paste("PG",pg_id,"l2fc mean\nSK=",format(sk,digits=2),"SD=",format(sd(cmean),digits=2)))

    temp <- fcmat2.se[sk.list,]
    cmean <- colMeans(temp)
    plot(density(cmean),xlim=c(0,0.4),xlab="mean(gene-wise se)",ylim=c(0,40),cex.lab=0.8,cex.axis=0.8,cex.main=0.8,
         main=paste("PG",pg_id,"SE mean"))

    temp <- fcmat2.pv[sk.list,]
    temp[temp<1e-3] <- 1e-3
    temp <- log10(temp)
    cmean <- colMeans(temp)
    cmean <- cmean[cmean<0]
    plot(density(cmean),xlim=c(-3,0),xlab="mean(gene-wise pv)",ylim=c(0,0.1),cex.lab=0.8,cex.axis=0.8,cex.main=0.8,
         main=paste("PG",pg_id,"PV mean"))

    do.median <- F
    if(do.median) {
      temp <- fcmat2[sk.list,]
      cmed <- apply(temp,FUN=median,MARGIN=2)
      plot(density(cmed),xlim=c(-0.5,0.5),xlab="median(gene-wise l2fc)",ylim=c(0,15),cex.lab=0.8,cex.axis=0.8,cex.main=0.8,
           main=paste("PG",pg_id,"l2fc median"))

      temp <- fcmat2.se[sk.list,]
      cmed <- apply(temp,FUN=median,MARGIN=2)
      plot(density(cmed),xlim=c(0,0.4),xlab="median(gene-wise se)",ylim=c(0,50),cex.lab=0.8,cex.axis=0.8,cex.main=0.8,
           main=paste("PG",pg_id,"SE median"))

      temp <- fcmat2.pv[sk.list,]
      cmed <- apply(temp,FUN=median,MARGIN=2)
      cmed[cmed<1e-6] <- 1e-6
      cmed <- log10(cmed)
      cmed <- cmed[cmed<0]
      plot(density(cmed),xlim=c(-6,0),xlab="median(gene-wise pv)",ylim=c(0,40),cex.lab=0.8,cex.axis=0.8,cex.main=0.8,
           main=paste("PG",pg_id,"PV median"))
    }
    if(!to.file) browser()
  }
  if(to.file) dev.off()

 }
