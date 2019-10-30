#--------------------------------------------------------------------------------------
#' Create a caalog of all of the possible pathways and wruite this to a
#' file that can be used to hand select specific pathways to use
#'
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
pathwayCatalog = function(){
  printCurrentFunction()
  load("../input/processed_pathway_data/msigdb_PATHWAYS.RData")
  load("../input/processed_pathway_data/RYAN_PATHWAYS.RData")
  load("../input/processed_pathway_data/bioplanet_PATHWAYS.RData")

  all_pathways <- rbind(bioplanet_PATHWAYS,msigdb_PATHWAYS,RYAN_PATHWAYS)
  all_pathways <- all_pathways[,1:4]
  all_pathways$superclass <- "-"
  all_pathways$useme <- 0

  all_pathways <- all_pathways[c("pathway","useme","superclass","pathset","pathway_class","path_info")]
  file <- "../input/processed_pathway_data/pathway_catalog.xlsx"
  write.xlsx(all_pathways,file)
}
