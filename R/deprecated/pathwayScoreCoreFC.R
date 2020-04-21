#' Pathway Score Core - FC
#'
#' Computes fold change pathway scores.
#'
#' This fast implementation of fold change pathway scores uses matrix
#' multiplication. The score is simply: mean(fold change of genes in pathway) -
#' mean(fold change of genes outside pathway).
#'
#' @param fcdata Sample by gene matrix of log2(fold change)'s. Rownames are
#'   sample keys and colnames are genes.
#' @param pathset Name of pathway set.
#' @param dataset Name of data set.
#' @param chem_dict Dataframe with one row per sample key and seven columns:
#'   sample_key, sample_id, conc, time, casrn, name, dtxsid.
#' @param pathway_data Named ist of gene name vectors. Each element is one
#'   pathway, defined by the genes it contains.
#'
#' @import openxlsx
#'
#' @return Dataframe with one row per chemical/conc/pathway combination. Columns
#'   are: sample_id, dtxsid, casrn, name, time, conc, pathset,
#'   pathway, size (pathway size accounting for missing values), mean_fc_scaled_in,
#'   mean_fc_scaled_out, pathway_score.
#' @export
pathwayScoreCoreFC <- function(fcdata,
                               pathset,
                               dataset,
                               chem_dict,
                               pathway_data) {

  #Strategy: Matrix multiplication against logical in/out matrix gives fc sum over pathway genes

  #inmat is a logical gene by pathway matrix with TRUE when a a gene is in the pathway
  # and false otherwise
  allgenes = colnames(fcdata)
  inmat = sapply(pathway_data, function(x){allgenes %in% x})

  #Change NAs to zero in data to avoid changing sum
  existingdata = fcdata
  existingdata[is.na(fcdata)] = 0
  #inlengths is a sample by pathway matrix that gives the number of non-missing genes present
  # in a pathway
  inlengths = (!is.na(fcdata)) %*% inmat
  insums = existingdata %*% inmat #insums is a sample by pathway matrix of l2fc sums in pathway
  inmeans = insums/inlengths  #Divide by length of actual number of genes used for mean

  #repeat for fold changes outside the pathway
  outsums = existingdata %*% (!inmat)
  outlengths = (!is.na(fcdata)) %*% (!inmat)
  outmeans = outsums/outlengths

  pathwayscore = inmeans - outmeans #final score is mean in - mean out

  #set up output columns
  name.list <- c("sample_id",
                 "dtxsid",
                 "casrn",
                 "name",
                 "time",
                 "conc",
                 "pathset",
                 "pathway",
                 "size",
                 "mean_fc_scaled_in",
                 "mean_fc_scaled_out",
                 "pathway_score")

  #pathscoremat will be arranged as every pathway for a given sample key,
  # then every pathway for the next sample, etc.
  pathscoremat = as.data.frame(matrix(0,nrow = nrow(pathwayscore)*ncol(pathwayscore),
                                      ncol = length(name.list)), stringsAsFactors = F)
  colnames(pathscoremat) = name.list
  # fill in columns for output
  pathscoremat$pathway_score = as.vector(t(pathwayscore))
  pathscoremat$mean_fc_scaled_in = as.vector(t(inmeans))
  pathscoremat$mean_fc_scaled_out = as.vector(t(outmeans))
  pathscoremat$size = as.vector(t(inlengths))
  pathscoremat$pathway = rep(names(pathway_data), nrow(fcdata))
  pathscoremat$pathset = rep(pathset, nrow(pathscoremat))

  #use chem_dict to attach chem information; assumes all id information is present in chem_dict
  sk.list = rownames(fcdata)
  chemkey = chem_dict[match(as.vector(sk.list), chem_dict$sample_key), name.list[1:6]]
  keyrepeat = rep(1:length(sk.list), each = length(pathway_data))
  pathscoremat[,1:6] = chemkey[keyrepeat,]

  return(pathscoremat)
  cat("Score Calculated\n")
}
