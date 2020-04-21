#--------------------------------------------------------------------------------------
#' Make heatmaps of the FCMAT and FCMAT.SE by plate groups
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
noiseHunter1.hm <- function(to.file=F,
                               do.load=F,
                               dataset="DMEM_6hr_screen_normal_pe_1",
                               sigset="pilot_large_all_CMAP",
                               sigcatalog="signatureDB_master_catalog 2020-03-10",
                               method="mygsea",
                               cutoff=0.25,
                            ngene=-1){
  printCurrentFunction()
  basedir="../input/fcdata/"
  if(do.load) {
    file <- paste("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData",sep="")
    print(file)
    load(file=file)
    signaturecr <<- SIGNATURE_CR
    annotations <<- signatureCatalogLoader(sigset,sigcatalog)

    file <- paste0("../input/chemicals/",dataset,"_chemical_map.xlsx")
    cmap <<- read.xlsx(file)

    file <- paste0(basedir,"FCMAT2_",dataset,".RData")
    print(file)
    load(file)
    FCMAT2 <<- FCMAT2
    file <- paste0(basedir,"FCMAT2.SE.",dataset,".RData")
    print(file)
    load(file)
    FCMAT2.SE <<- FCMAT2.SE

    file <- "../input/signatures/httr_mcf7_ph1_DMSO_by_gene.xlsx"
    genecounts <<- read.xlsx(file)
  }
  if(to.file) {
    fname <- paste0("../output/noiseHunter/noiseHunter1 ",cutoff," heatmap.pdf")
    if(ngene>0) fname <- paste0("../output/noiseHunter/noiseHunter1 ",cutoff," heatmap ngene ",ngene,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }

  fcmat2 <- FCMAT2
  fcmat2[is.na(fcmat2)] <- 0
  fcmat2.se <- FCMAT2.SE
  fcmat2.se[is.na(fcmat2.se)] <- 0
  gene.list <- colnames(fcmat2)

  gc <- unique(genecounts[,c("gene_symbol","log2_CPM_mean")])
  names(gc) <- c("gene","l2cpm")
  gc <- gc[!duplicated(gc$gene),]
  rownames(gc) <- gc$gene
  gc <- gc[gene.list,]
  x <- gc$l2cpm
  x[is.na(x)] <- 0
  gc$l2cpm <- x

  gc$color <- "white"
  gc[x>0,"color"] <- "black"
  gc[x>2,"color"] <- "gray"
  gc[x>5,"color"] <- "green"
  gc[x>10,"color"] <- "yellow"
  gc[x>15,"color"] <- "red"

  rownames(cmap) <- cmap$sample_key
  pg.list <- sort(unique(cmap$pg_id))
  cutoff.l2fcfc <- 2
  cutoff.se <- 0.75
  if(ngene== -1) ngene <- ncol(FCMAT2)
  for(i in 1:length(pg.list)) {
    pg_id <- pg.list[i]
    cat("Plate group",pg_id,"\n")
    sk.list <- cmap[cmap$pg_id==pg_id,"sample_key"]
    col.list <- gc$color

    temp <- fcmat2[sk.list,]
    temp <- temp[,1:ngene]
    col.list <- col.list[1:ngene]
    temp[temp>cutoff.l2fcfc] <- cutoff.l2fcfc
    temp[temp< -cutoff.l2fcfc] <- -cutoff.l2fcfc
    result <- heatmap.2(temp,
                        margins=c(10,10),
                        dendrogram="col",
                        scale="none",
                        main=paste("FCMAT2 plate group:",pg_id),
                        xlab="",
                        ylab="",
                        cexCol=0.1,
                        cexRow=0.1,
                        Rowv=F,
                        Colv=T,
                        trace="none",
                        hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                        key=T,
                        col=brewer.pal(9,"PiYG"),
                        key.title="Key",
                        key.xlab="l2fc",
                        cex.main=1,
                        ColSideColors=col.list)

    if(!to.file) browser()

    temp <- fcmat2.se[sk.list,]
     temp <- temp[,1:ngene]
    col.list <- col.list[1:ngene]
    temp[temp>cutoff.se] <- cutoff.se
    temp[temp< -cutoff.se] <- -cutoff.se
    result <- heatmap.2(temp,
                        margins=c(5,5),
                        dendrogram="col",
                        scale="none",
                        main=paste("FCMAT2 SE plate group:",pg_id),
                        xlab="",
                        ylab="",
                        cexCol=0.1,
                        cexRow=0.1,
                        Rowv=F,
                        Colv=T,
                        trace="none",
                        hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                        key=T,
                        col=brewer.pal(9,"Reds"),
                        key.title="Key",
                        key.xlab="l2fc",
                        cex.main=1,
                        ColSideColors=col.list)

    if(!to.file) browser()
  }


  if(to.file) dev.off()
}
