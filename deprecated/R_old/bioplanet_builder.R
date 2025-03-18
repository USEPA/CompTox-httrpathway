#' BioPlanet Builder
#'
#' Converts BioPlanet data into usable pathway data.
#'
#' This function shows how BioPlanet data was converted to usable pathway files.
#' As BioPlanet is updated, this function will have to be updated. It requires
#' two downloaded .csv files with location specified by pathfile and catfile.
#' It saves usable pathway files with location specified by pwayout and
#' pdataout to disk.
#'
#' @param pathfile File name of bioplanet_pathway.csv.
#' @param catfile File name of bioplanet_pathway_category.csv.
#' @param pwayout File name of bioplanet_PATHWAYS.RDS
#' @param pdataout File name to save to
#'
#' @return No output.
#' @export bioplanet_builder
#'
#' @importFrom stats aggregate
#' @importFrom utils read.csv

bioplanet_builder = function(pathfile = "../input/processed_pathway_data/bioplanet_pathway.csv",
                             catfile = "../input/processed_pathway_data/bioplanet_pathway_category.csv",
                             pwayout = "../input/processed_pathway_data/bioplanet_PATHWAYS.RDS",
                             pdataout = "../input/processed_pathway_data/PATHWAY_LIST_bioplanet.RDS"){

  pways = read.csv(pathfile, stringsAsFactors = F)
  pwaycats = read.csv(catfile, stringsAsFactors = F)

  #an ad hoc solution to the bizarre fact that bioplanet has two different pathways with the same name
  pways$PATHWAY_NAME[pways$PATHWAY_ID == "bioplanet_703"] = "p75 neurotrophin receptor-mediated signaling_manual"
  pwaycats$PATHWAY_NAME[pwaycats$PATHWAY_ID == "bioplanet_703"] = "p75 neurotrophin receptor-mediated signaling_manual"

  ls = aggregate(GENE_SYMBOL ~ PATHWAY_NAME,pways, length)
  pways = aggregate(GENE_SYMBOL ~ PATHWAY_NAME,pways, function(x){paste0(x, collapse = "|")})
  pways$ngene = ls$GENE_SYMBOL
  pwaycats = pwaycats[pwaycats$PATHWAY_NAME != "",]
  pwaycats = aggregate(CATEGORY_NAME ~ PATHWAY_NAME,pwaycats, function(x){paste0(x, collapse = "|")})

  pways = cbind(pways,pwaycats$CATEGORY_NAME[match(pways$PATHWAY_NAME, pwaycats$PATHWAY_NAME)])

  #this naming part of the code is very brittle, so exercise caution if inputs or previous code are changed
  colnames(pways) = c("pathway", "gene_list", "ngene","path_info")
  pways$pathset = "bioplanet"
  pways$pathway_class = NA

  bioplanet_PATHWAYS = pways[,c(1,5,6,4,2,3)]
  #save(bioplanet_PATHWAYS, file = pwayout)
  saveRDS(bioplanet_PATHWAYS, pwayout)

  pathway_data = strsplit(bioplanet_PATHWAYS$gene_list, "\\|")
  names(pathway_data) = bioplanet_PATHWAYS$pathway
  #save(pathway_data, file = pdataout)
  saveRDS(pathway_data, pdataout)
}
