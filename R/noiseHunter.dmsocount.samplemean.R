#source("https://bioconductor.org/biocLite.R")
#.libPaths( "C:/Program Files/R/R-3.3.3/library")
library(fastmatch)
library(reshape2)
#source("R/general_utils.R")
#options(java.parameters = "-Xmx8000m")
#library(httRws)
library(openxlsx)
library(tidyr)
library(grDevices)
library(RColorBrewer)
library(gplots)
library(ape)
library(pracma)
#--------------------------------------------------------------------------------------
#' Find issues with the DMSO samples that may be causing the odd noise issue
#'
#' @param dir The directory to work in
#'
#--------------------------------------------------------------------------------------
noiseHunter.dmsocount.samplemean <- function(to.file=F,
                                             do.read=F,
                                             dirin="../input/httr_mcf7_screen/raw_ctrl/",
                                             welltype="dmso",
                                             well_type_prefix="",
                                             dataset="DMEM_6hr_screen_normal_pe_1") {
  printCurrentFunction()
  if(do.read) {
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

    ds2 <- str_replace(dataset,"_pgnorm","")
    file <- paste0("../input/chemicals/",ds2,"_chemical_map.xlsx")
    cmap <<- read.xlsx(file)
    file <- paste0("../input/fcdata/FCMAT2_",dataset,".RData")
    print(file)
    load(file)
    FCMAT2 <<- FCMAT2
  }
  if(to.file) {
    fname <- paste0("../output/noiseHunter/noiseHunter1 ",dataset," dmsocount.samplemean.pdf")
    pdf(file=fname,width=7,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(4,3),mar=c(4,4,2,2))
  dmso.all <- CTRLMAT
  genes <- separate(dmso.all,"probe_id",sep="_",into=c("genes","id"))
  dmso.all$gene <- genes[,"genes"]

  fcmat2 <- FCMAT2
  fcmat2[is.na(fcmat2)] <- 0
  gene.list <- colnames(fcmat2)
  gene.list <- colnames(fcmat2)
  gene.list <- gene.list[is.element(gene.list,dmso.all$gene)]
  fcmat2 <- fcmat2[,gene.list]

  rownames(cmap) <- cmap$sample_key
  pg.list <- sort(unique(cmap$pg_id))

  for(i in 1:length(pg.list)) {
    pg_id <- pg.list[i]
    cat("Plate group",pg_id,"\n")
    sk.list <- cmap[cmap$pg_id==pg_id,"sample_key"]
    temp <- fcmat2[sk.list,]

    dmso <- dmso.all[is.element(dmso.all$pg_id,pg_id),]
    res <- reshape2::dcast(dmso,sample_id~gene,value.var="probe_count",fill=0.5,fun.aggregate=max)
    rownames(res) <- res[,1]
    res <- res[,2:ncol(res)]

    glist <- colnames(temp)
    glist <- glist[is.element(glist,names(res))]
    cat("   genes:",length(glist),"\n")
    res <- res[,glist]
    temp <- temp[,glist]
    cmean <- colMeans(temp)
    cmed <- apply(temp,FUN=median,MARGIN=2)

    rs <- rowSums(res)
    for(i in 1:nrow(res)) res[i,] <- res[i,]*1000000/rs[i]
    res <- log2(res)
    dmean <- colMeans(res)
    dmed <- apply(res,FUN=median,MARGIN=2)
    sp <- cor.test(dmean,cmean,method="spearman",exact=F)
    rho <- sp$estimate
    pe <- cor.test(dmean,cmean,method="pearson",exact=F)
    r2 <- pe$estimate



    main <- paste("PG:",pg_id,":",format(rho,digits=2),":",format(r2,digits=2))
    plot(cmean~dmean,xlab="Baseline log2(CPM)",ylab="l2fc",main=main,
         cex.lab=1.2,cex.axis=1.2,pch=".",xlim=c(0,12),ylim=c(-0.4,0.4),col="gray",bg="gray")
    lines(c(0,100),c(0,0),col="red")
    lines(c(0,0),c(-1,1),col="red")
    x <- density(dmean)
    lines(x$y~x$x,col="red")

    x <- density(cmean)
    yval <- x$x
    xval <- x$y
    lines(yval~xval,col="red")

    mod <- lm(cmean~dmean)
    slope <- mod$coefficients[2]
    intercept <- mod$coefficients[1]
    x1 <- 0
    x2 <- 12
    y1 <- intercept
    y2 <- intercept + slope*(x2-x1)
    lines(c(x1,x2),c(y1,y2),col="blue")
    if(!to.file) browser()
  }

  if(to.file) dev.off()
}
