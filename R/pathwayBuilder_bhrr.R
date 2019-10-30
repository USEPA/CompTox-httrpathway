#' bhrr Pathway Builder
#'
#' Builds bhrr pathset based on msigdb, bioplanet, and ryan pathsets.
#'
#' Shows how the bhrr pathset was built from pre-existing pathsets.
#'
#' @return No output.
#' @export
pathwayBuilder_bhrr = function(){
  load("../input/processed_pathway_data/msigdb_PATHWAYS.RData")

  hallmark = msigdb_PATHWAYS[msigdb_PATHWAYS$pathway_class == "H: hallmark gene sets",]
  response = msigdb_PATHWAYS[msigdb_PATHWAYS$pathway_class == "C2: chemical and genetic perturbations",]
  response = response[grepl("RESPONSE",response$pathway),]

  load("../input/processed_pathway_data/RYAN_PATHWAYS.RData")
  load("../input/processed_pathway_data/bioplanet_PATHWAYS.RData")

  bhrr_PATHWAYS = rbind(bioplanet_PATHWAYS, hallmark, response, RYAN_PATHWAYS)
  save(bhrr_PATHWAYS, file = "input/processed_pathway_data/bhrr_PATHWAYS.RData")

  pathway_data = strsplit(bhrr_PATHWAYS$gene_list, "\\|")
  names(pathway_data) = bhrr_PATHWAYS$pathway
  save(pathway_data, file = "../input/processed_pathway_data/PATHWAY_LIST_bhrr.RData")
}
