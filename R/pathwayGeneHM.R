#--------------------------------------------------------------------------------------
#'
#' Build lane plots by chemical list and pathway class, across the datasets
#' @param to.file If TRUE, write plots to a file
#--------------------------------------------------------------------------------------
pathwayGeneHM <- function(to.file=F,
                          dataset="DMEM_6hr_pilot_normal_pe_1",
                          chemical.target="ER",
                          pathway.super_class="estrogen",
                          pathset="PathwaySet_20191107",
                          method = "gsva",
                          threshold=0.5) {
  printCurrentFunction()
  file <- "../input/chemicals/HTTR Pilot chemical annotation.xlsx"
  chems0 <- read.xlsx(file)
  chems0 <- chems0[is.element(chems0$target_key,chemical.target),]
  dtxsid.list <- chems0$dtxsid
  nchem <- length(dtxsid.list)

  file <- "../input/processed_pathway_data/pathway_catalog 2019-11-07.xlsx"
  catalog <- read.xlsx(file)
  catalog <- catalog[catalog$useme==1,]
  catalog <- catalog[is.element(catalog$super_class,pathway.super_class),]
  pathway.list <- catalog$pathway
  pathway.list <- sort(pathway.list)
  #pathway.list <- pathway.list[1:3]
  if(to.file) {
    fname <- paste0("../output/pathway_gene_heatmaps/pathwayGeneHM_",dataset,"_",chemical.target,"_",pathway.super_class,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }

  file <- paste0("../input/fcdata/FCMAT2_",dataset,".RData")
  print(file)
  load(file=file)
  fcmat2 <- FCMAT2

  file <- paste0("../input/fcdata/FCMAT2.PV.",dataset,".RData")
  print(file)
  load(file=file)
  fcmat2.pv <- FCMAT2.PV

  file <- paste0("../input/fcdata/FCMAT2.FCoverSE.",dataset,".RData")
  print(file)
  load(file=file)
  fcmat2.fcoverse <- FCMAT2.FCoverSE

  file <- paste0("../input/fcdata/CHEM_DICT_",dataset,".RData")
  load(file=file)
  chems <- CHEM_DICT
  chems <- chems[is.element(chems$dtxsid,dtxsid.list),]
  chemnames <- paste(chems$name,chems$conc)
  file <- paste0("../input/processed_pathway_data/PATHWAY_CATALOG_",pathset,".RData")
  load(file=file)
  file <- paste0("../input/processed_pathway_data/PATHWAY_GENE_LIST_",pathset,".RData")
  load(file=file)


  for(pathway in pathway.list) {
    cat(pathway,"\n")

    gene.list <- pathway_data[pathway][[1]]
    gene.list <- unique(gene.list)
    gene.list <- gene.list[is.element(gene.list,colnames(fcmat2))]

    for(mode in c("fc","fcse","pv")) {
      if(mode=="fc") mat <- fcmat2
      else if(mode=="pv") mat <- -log10(fcmat2.pv)
      else if(mode=="fcse") mat <- fcmat2.fcoverse

      mat2 <- mat[,gene.list]
      mat2[is.na(mat2)] <- 0

      mat3 <- mat2[chems$sample_key,]
      print(dim(mat3))
      if(mode=="fc") cutoff <- 2
      if(mode=="fcse") cutoff <- 5
      if(mode=="pv") cutoff <- 4
      mat3[mat3>cutoff] <- cutoff
      mat3[mat3< -cutoff] <- -cutoff
      rowsep <- NULL
      lastname <- ""
      for(i in 1:nrow(chems)) {
        if(chems[i,"name"]!=lastname) {
          lastname <- chems[i,"name"]
          rowsep <- c(rowsep,i)
        }
      }
      rowsep <- rowsep-1
      if(is.element(mode,c("fc","fcse"))) colset <- brewer.pal(9,"RdBu")
      else colset <- brewer.pal(9,"Reds")
      result <- heatmap.2(as.matrix(mat3),
                          margins=c(5,5),
                          scale="none",
                          main=paste(mode,"\n",pathway),
                          xlab="",
                          ylab="",
                          cexCol=0.1,
                          cexRow=0.3,
                          col=colset,
                          Rowv=F,
                          Colv=T,
                          dendrogram="column",
                          trace="none",
                          hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                          key=T,
                          key.title="Key",
                          key.xlab="l2fc",
                          labRow=chemnames,
                          rowsep=rowsep,
                          sepcolor="gray",
                          cex.main=1)

      if(!to.file) browser()
    }
  }
  if(to.file) dev.off()
}

