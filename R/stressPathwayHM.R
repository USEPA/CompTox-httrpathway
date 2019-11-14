#--------------------------------------------------------------------------------------
#'
#' Build a heatmap of the stress genes
#' @param to.file If TRUE, write plots to a file
#--------------------------------------------------------------------------------------
stressPathwayHM <- function(to.file=F,
                            dataset="DMEM_6hr_pilot_normal_pe_0",
                            pathset="PathwaySet_20191107",
                            threshold=0.5,
                            method="fc") {
  printCurrentFunction()
  if(to.file) {
    fname <- paste0("../output/miscplots/stressPathwayHM ",dataset,"_",pathset,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }

  file <- paste0("../input/fcdata/FCMAT2_",dataset,".RData")
  load(file=file)
  mat <- FCMAT2
  file <- paste0("../input/fcdata/CHEM_DICT_",dataset,".RData")
  load(file=file)
  chems <- CHEM_DICT
  chemnames <- paste(chems$name,chems$conc)
  file <- paste0("../input/processed_pathway_data/PATHWAY_CATALOG_",pathset,".RData")
  load(file=file)
  file <- paste0("../input/processed_pathway_data/PATHWAY_GENE_LIST_",pathset,".RData")
  load(file=file)

  sclass.list <- c(
    "dna damage" ,
    "cell cycle" ,
    "apoptosis" ,
    "hdac",
    "nfkb",
    "tgfb",
    "egfr",
    "tnf" ,
    "fgfr" ,
    "cytoskeleton",
    "oxidative stress" ,
    "cytotoxicity" ,
    "mitochondria",
    "proliferation",
    "microtubule",
    "hypoxia",
    "heat shock",
    "hif1a",
    "h1f2a",
    "sterol processing",
    "estrogen",
    "androgen",
    "adipogenesis",
    "fatty acid",
    "p450",
    "Amiodarone",
    "conazole",
    "Cyproterone",
    "fibrate",
    "Flutamide",
    "glitazone",
    "Nilutamide",
    "Rotenone",
    "statin",
    "Tamoxifen",
    "Testosterone",
    "Trichostatin A"
    )

  #sclass.list <- sclass.list[1:3]
  for(sclass in sclass.list) {
    pathway.list <- pathway_catalog[is.element(pathway_catalog$super_class,sclass),"pathway"]
    gene.list <- NULL
    for(pathway in pathway.list) gene.list <- c(gene.list,pathway_data[pathway][[1]])
    gene.list <- unique(gene.list)
    gene.list <- gene.list[is.element(gene.list,colnames(mat))]
    cat(sclass,": ",length(gene.list),"\n")
    mat2 <- mat[,gene.list]
    mat2[is.na(mat2)] <- 0
    mat2[abs(mat2)<threshold] <- 0
    cs <- colSums(mat2)
    mat2 <- mat2[,cs>0]

    cutoff <- 2
    mat[mat>cutoff] <- cutoff
    mat[mat< -cutoff] <- -cutoff
    cat(sclass,"\n")
    print(dim(mat2))
    rowsep <- seq(from=0,to=,44*8,by=8)

    rowsep <- NULL
    lastname <- ""
    for(i in 1:nrow(chems)) {
      if(chems[i,"name"]!=lastname) {
        lastname <- chems[i,"name"]
        rowsep <- c(rowsep,i)
      }
    }
    rowsep <- rowsep-1

    result <- heatmap.2(as.matrix(mat2),
                        margins=c(5,5),
                        scale="none",
                        main=sclass,
                        xlab="",
                        ylab="",
                        cexCol=0.1,
                        cexRow=0.1,
                        col=brewer.pal(9,"RdBu"),
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
                        sepcolor="gray")

    if(!to.file) browser()
  }
  if(to.file) dev.off()
}

