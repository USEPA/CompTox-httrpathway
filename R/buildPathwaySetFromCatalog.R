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
buildPathwaySetFromCatalog = function(catalog.file="../input/processed_pathway_data/pathway_catalog 2019-10-31.xlsx",
                                      pathsetname="PathwaySet_20191031",do.filter=T){
  printCurrentFunction()
  load("../input/processed_pathway_data/msigdb_PATHWAYS.RData")
  load("../input/processed_pathway_data/RYAN_PATHWAYS.RData")
  load("../input/processed_pathway_data/bioplanet_PATHWAYS.RData")

  all_pathways <- rbind(bioplanet_PATHWAYS,msigdb_PATHWAYS,RYAN_PATHWAYS)
  rownames(all_pathways) <- all_pathways$pathway
  catalog <- read.xlsx(catalog.file)
  if(do.filter) catalog <- catalog[catalog$useme>0,]
  pathway.list <- catalog$pathway
  pathways <- all_pathways[pathway.list,]
  pathways$superclass <- catalog$superclass

  save(pathways, file = paste0("../input/processed_pathway_data/PATHWAYS_",pathsetname,".RData"))

  pathway_data = strsplit(pathways$gene_list, "\\|")
  print(length(pathway_data))
  names(pathway_data) = pathways$pathway
  save(pathway_data, file = paste0("../input/processed_pathway_data/PATHWAY_LIST_",pathsetname,".RData"))
}
