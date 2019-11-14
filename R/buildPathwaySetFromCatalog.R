#--------------------------------------------------------------------------------------
#' Create a caalog of all of the possible pathways and wruite this to a
#' file that can be used to hand select specific pathways to use
#'
#' @param catalog.file This is the name of the catalog file
#' @param pathsetname THis is the name of the pathway set and will be used for the
#' name of the output file
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
buildPathwaySetFromCatalog = function(catalog.file="../input/processed_pathway_data/pathway_catalog 2019-11-07.xlsx",
                                      pathset="PathwaySet_20191107",
                                      dataset="DMEM_6hr_pilot_normal_pe_0",
                                      do.filter=T,
                                      nrandom=500,
                                      mean.random=100,
                                      sd.random=40){
  printCurrentFunction()
  load("../input/processed_pathway_data/msigdb_PATHWAYS.RData")
  load("../input/processed_pathway_data/RYAN_PATHWAYS.RData")
  load("../input/processed_pathway_data/bioplanet_PATHWAYS.RData")
  load("../input/processed_pathway_data/CMAP_PATHWAYS.RData")

  all_pathways <- rbind(bioplanet_PATHWAYS,msigdb_PATHWAYS,RYAN_PATHWAYS,CMAP_PATHWAYS)
  rownames(all_pathways) <- all_pathways$pathway
  catalog <- read.xlsx(catalog.file)
  rownames(catalog) <- catalog$pathway
  if(do.filter) catalog <- catalog[catalog$useme>0,]
  pathway.list <- catalog$pathway
  pathways <- all_pathways[pathway.list,]
  pathways$super_class <- catalog$super_class

  load(file=paste0("../input/fcdata/FCMAT2_",dataset,".RData"))
  gene.list <- colnames(FCMAT2)

  #for(i in 1:nrow(pathways)) {
  #  pathway <- pathways[i,"pathway"]
  #  glist <- all_pathways[pathway,"gene_list"]
  #  glist <- str_split(glist,"\\|")[[1]]
  #  ng <- length(glist[is.element(glist,gene.list)])
  #  pathways[i,"ngene_in_data"] <- ng
  #  pathways[i,"pclass"] <- catalog[pathway,"pclass"]
  #}

  if(nrandom>0) {
    name.list <- names(pathways)
    prandom <- as.data.frame(matrix(nrow=nrandom,ncol=length(name.list)))
    names(prandom) <- name.list
    n <- 0
    for(i in 1:nrandom) {
      while(n<10) n <- round(rnorm(1,mean=mean.random,sd=sd.random))

      glist <- sample(gene.list,n)
      pname <- paste0("Random_",i)
      prandom[i,"pathway"] <- pname
      prandom[i,"pathset"] <- "random"
      prandom[i,"pathway_class"] <- "random"
      prandom[i,"path_info"] <- "random"
      prandom[i,"gene_list"] <- paste(glist,collapse="|")
      prandom[i,"ngene"] <- n
      prandom[i,"super_class"] <- "random"
    }
    pathways <- rbind(pathways,prandom)
  }

  pathway_data = strsplit(pathways$gene_list, "\\|")
  print(length(pathway_data))
  names(pathway_data) = pathways$pathway
  save(pathway_data, file = paste0("../input/processed_pathway_data/PATHWAY_GENE_LIST_",pathset,".RData"))

  name.list <- names(pathways)
  name.list <- name.list[!is.element(name.list,"gene_list")]
  pathway_catalog <- pathways[,name.list]
  save(pathway_catalog, file = paste0("../input/processed_pathway_data/PATHWAY_CATALOG_",pathset,".RData"))
}
