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
#' @param pwayout File name of bioplanet_PATHWAYS.RData
#' @param pdataout File name of
#'
#' @return No output.
#' @export
#'
#' @import stats
#' @import utils
bioplanet_builder = function(pathfile = "../input/processed_pathway_data/bioplanet_pathway.csv",
                             catfile = "../input/processed_pathway_data/bioplanet_pathway_category.csv",
                             pwayout = "../input/processed_pathway_data/bioplanet_PATHWAYS.RData",
                             pdataout = "../input/processed_pathway_data/PATHWAY_LIST_bioplanet.RData"){

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
  save(bioplanet_PATHWAYS, file = pwayout)

  pathway_data = strsplit(bioplanet_PATHWAYS$gene_list, "\\|")
  names(pathway_data) = bioplanet_PATHWAYS$pathway
  save(pathway_data, file = pdataout)
}
