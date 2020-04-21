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
noiseHunter1.pgplots.filter <- function(to.file=F,
                                 do.load=F,do.load.ctrl=F,
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
  if(do.load.ctrl) {
    dirin <- "../input/httr_mcf7_screen/raw_ctrl/"
    welltype <- "dmso"
    well_type_prefix <- ""
    file <- paste0(dirin,"httr-ph-i-raw-pl-bl_",welltype,"-v1.tsv")
    print(file)
    mat <- read.delim(file,stringsAsFactors=F)
    x <- unite(mat[,c("plate_id","block_id","pg_id","well_id")],"rowkey")
    mat$rowkey <- x[,1]
    CTRLMAT <<- mat
    pmat <- unique(CTRLMAT[,c("plate_id","block_id","pg_id","well_id")])
    x <- unite(pmat,"rowkey")
    rownames(pmat) <- x[,1]
    PLATEMAT <<- pmat
  }
  dmso.all <- CTRLMAT
  genes <- separate(dmso.all,"probe_id",sep="_",into=c("genes","id"))
  dmso.all$gene <- genes[,"genes"]

  #---------------------------------------------------------------------------------
  # se
  #---------------------------------------------------------------------------------
  if(to.file) {
    fname <- paste0("../output/noiseHunter/noiseHunter ",dataset," pg distribution filtered.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(6,4),mar=c(4,3,3,2))

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
    dmso <- dmso.all[is.element(dmso.all$pg_id,pg_id),]
    res <- reshape2::dcast(dmso,sample_id~gene,value.var="probe_count",fill=0.5,fun.aggregate=max)
    rownames(res) <- res[,1]
    res <- res[,2:ncol(res)]
    res <- as.matrix(res)
    gene.list <- colnames(FCMAT2)
    gene.list <- gene.list[is.element(gene.list,colnames(res))]
    res2 <- res[,gene.list]
    rs <- rowSums(res2)
    for(i in 1:nrow(res2)) res2 <- res2*1000000/rs[i]
    res2 <- log2(res2)
    cmdmso <- colMeans(res2)
    cat("Plate group",pg_id,"\n")
    sk.list <- cmap[cmap$pg_id==pg_id,"sample_key"]

    temp <- fcmat2[sk.list,gene.list]
    temp0 <- temp
    stemp <- fcmat2.se[sk.list,gene.list]
    for(cutoff in c(0.8,0.6,0.4,0.2)) {
      mask <- stemp
      mask[] <- 1
      mask[stemp>cutoff] <- 0
      count <- sum(as.numeric(as.matrix(mask)))
      ratio <- count / nrow(mask) / ncol(mask)
      ratio <- 1-ratio
      temp[mask==0] <- NA
      cmean <- colMeans(temp,na.rm=T)
      plot(density(cmean),xlim=c(-1,1),xlab="mean(gene-wise l2fc)",ylim=c(0,15),cex.lab=0.8,cex.axis=0.8,cex.main=0.8,
           main=paste("PG",pg_id,"l2fc mean se<",cutoff,"\n frac removed:",format(ratio,digits=2)))
    }
    for(cutoff in c(0,2,4,5)){
      temp <- temp0
      cmdmso
      temp <- temp[,cmdmso>cutoff]
      #print(dim(temp))
      doit <- T
      if(is.null(dim(temp))) doit <- F
      else if(ncol(temp)==0) doit <- F
      #if(pg_id==26) browser()
      if(doit) {
        cmean <- colMeans(temp,na.rm=T)
        plot(density(cmean),xlim=c(-1,1),xlab="mean(gene-wise l2fc)",ylim=c(0,15),cex.lab=0.8,cex.axis=0.8,cex.main=0.8,
             main=paste("PG",pg_id,"l2fc mean |log2 DMSO|>",cutoff))
      }
      else {
        plot(c(0,0),main="missing data")
      }
    }
    if(!to.file) browser()
  }
  if(to.file) dev.off()

}
