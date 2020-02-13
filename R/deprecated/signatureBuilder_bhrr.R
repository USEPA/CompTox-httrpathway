#' bhrr Pathway Builder
#'
#' Builds bhrr pathset based on msigdb, bioplanet, and ryan pathsets.
#'
#' Shows how the bhrr pathset was built from pre-existing pathsets.
#'
#' @return No output.
#' @export
signatureBuilder_bhrr = function(){
  load("../input/processed_signature_data/msigdb_PATHWAYS.RData")

  hallmark = msigdb_PATHWAYS[msigdb_PATHWAYS$signature_class == "H: hallmark gene sets",]
  response = msigdb_PATHWAYS[msigdb_PATHWAYS$signature_class == "C2: chemical and genetic perturbations",]
  response = response[grepl("RESPONSE",response$signature),]

  load("../input/processed_signature_data/RYAN_PATHWAYS.RData")
  load("../input/processed_signature_data/bioplanet_PATHWAYS.RData")

  bhrr_PATHWAYS = rbind(bioplanet_PATHWAYS, hallmark, response, RYAN_PATHWAYS)
  save(bhrr_PATHWAYS, file = "input/processed_signature_data/bhrr_PATHWAYS.RData")

  signature_data = strsplit(bhrr_PATHWAYS$gene_list, "\\|")
  names(signature_data) = bhrr_PATHWAYS$signature
  save(signature_data, file = "../input/processed_signature_data/PATHWAY_LIST_bhrr.RData")
}
