#' Signature Score Core - FC
#'
#' Computes fold change signature scores.
#'
#' This fast implementation of fold change signature scores uses matrix
#' multiplication. The score is simply: mean(fold change of genes in signature) -
#' mean(fold change of genes outside signature).
#'
#' @param fcdata Sample by gene matrix of log2(fold change)'s. Rownames are
#'   sample keys and colnames are genes.
#' @param sigset Name of signature set.
#' @param dataset Name of data set.
#' @param chem_dict Dataframe with one row per sample key and seven columns:
#'   sample_key, sample_id, conc, time, casrn, name, dtxsid.
#' @param signature_data Named ist of gene name vectors. Each element is one
#'   signature, defined by the genes it contains.
#'
#' @import openxlsx
#'
#' @return Dataframe with one row per chemical/conc/signature combination. Columns
#'   are: sample_id, dtxsid, casrn, name, time, conc, sigset,
#'   signature, size (signature size accounting for missing values), mean_fc_scaled_in,
#'   mean_fc_scaled_out, signature_score.
#' @export
signatureScoreCoreFC <- function(fcdata,
                               sigset,
                               dataset,
                               chem_dict,
                               signature_data) {

  #Strategy: Matrix multiplication against logical in/out matrix gives fc sum over signature genes

  #inmat is a logical gene by signature matrix with TRUE when a a gene is in the signature
  # and false otherwise
  allgenes = colnames(fcdata)
  inmat = sapply(signature_data, function(x){allgenes %in% x})

  #Change NAs to zero in data to avoid changing sum
  existingdata = fcdata
  existingdata[is.na(fcdata)] = 0
  #inlengths is a sample by signature matrix that gives the number of non-missing genes present
  # in a signature
  inlengths = (!is.na(fcdata)) %*% inmat
  insums = existingdata %*% inmat #insums is a sample by signature matrix of l2fc sums in signature
  inmeans = insums/inlengths  #Divide by length of actual number of genes used for mean

  #repeat for fold changes outside the signature
  outsums = existingdata %*% (!inmat)
  outlengths = (!is.na(fcdata)) %*% (!inmat)
  outmeans = outsums/outlengths

  signaturescore = inmeans - outmeans #final score is mean in - mean out

  #set up output columns
  name.list <- c("sample_id",
                 "dtxsid",
                 "casrn",
                 "name",
                 "time",
                 "conc",
                 "sigset",
                 "signature",
                 "size",
                 "mean_fc_scaled_in",
                 "mean_fc_scaled_out",
                 "signature_score")

  #signaturescoremat will be arranged as every signature for a given sample key,
  # then every signature for the next sample, etc.
  signaturescoremat = as.data.frame(matrix(0,nrow = nrow(signaturescore)*ncol(signaturescore),
                                      ncol = length(name.list)), stringsAsFactors = F)
  colnames(signaturescoremat) = name.list
  # fill in columns for output
  signaturescoremat$signature_score = as.vector(t(signaturescore))
  signaturescoremat$mean_fc_scaled_in = as.vector(t(inmeans))
  signaturescoremat$mean_fc_scaled_out = as.vector(t(outmeans))
  signaturescoremat$size = as.vector(t(inlengths))
  signaturescoremat$signature = rep(names(signature_data), nrow(fcdata))
  signaturescoremat$sigset = rep(sigset, nrow(signaturescoremat))

  #use chem_dict to attach chem information; assumes all id information is present in chem_dict
  sk.list = rownames(fcdata)
  chemkey = chem_dict[match(as.vector(sk.list), chem_dict$sample_key), name.list[1:6]]
  keyrepeat = rep(1:length(sk.list), each = length(signature_data))
  signaturescoremat[,1:6] = chemkey[keyrepeat,]

  return(signaturescoremat)
  cat("Score Calculated\n")
}
