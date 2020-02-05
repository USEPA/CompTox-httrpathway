#' My Gene Set Enrichment Analysis
#' 
#' Performs tweaked version of single sample GSEA.
#' 
#' Based on the GSVA ssGSEA code. Main changes are: NAs are now handled correctly
#' and rank is now centered on zero instead of beginning at one. Since signature
#' sizes are undercounted here due to missing values, they are
#' assessed more accurately in signatureScoreCoreMYGSEA and limits are enforced
#' after scoring.
#'
#' @param X Transposed FCMAT2; i.e a gene by sample matrix of l2fc's including 
#'   rownames and colnames. Equivalent to expr in gsva.
#' @param geneSets Named list of signature definitions. Each element is a vector
#'   of gene names. Each element name is a signature name.E quivalent to 
#'   gset.idx.list in gsva.
#' @param min.sz Minimum signature size (deprecated).
#' @param max.sz Maximum signature size (deprecated)
#' @param alpha Power of R to use. Higher alpha will upweight more extreme
#'   ranks relative to middle ranks.
#' @param verbose verbose = T prints gene set length message.
#' @param useranks useranks = T uses ranks as in ssGSEA, while useranks = F
#'   uses the bare fold changes instead.
#' 
#' @importFrom GSVA filterGeneSets
#'
#' @return Outputs signature by sample matrix of signature scores.
#' @export
#'
#' @examples
#' geneSets = list(signature1 = c("ABC", "DEF"), signature2 = c("ABC", "GHI"))
#' X = matrix(c(1:3,3:1), nrow = 3)
#' colnames(X) = c("Sample1", "Sample2")
#' rownames(X) = c("ABC", "DEF", "GHI")
#' MYGSEA(X,geneSets)
MYGSEA = function(X, geneSets, min.sz = 1, max.sz = Inf, alpha = .25, verbose = T, useranks = T) {
  if (nrow(X) < 2) stop("Less than two genes in the input expression data matrix\n")
  
  mapped.geneSets <- lapply(geneSets, function(x, y) na.omit(match(x, y)), rownames(X))
  
  if (length(unlist(mapped.geneSets, use.names = FALSE)) == 0) {
    stop("No identifiers in the gene sets could be matched to the identifiers in the expression data.")
  }
  
  geneSets <- filterGeneSets(mapped.geneSets, min.sz = max(1, min.sz), max.sz = max.sz)
  if (length(geneSets) == 0) {
    stop("The gene set list is empty!  Filter may be too stringent.")
  }
  if (verbose) cat("Estimating MYGSEA scores for", length(geneSets), "gene sets.\n")

  p <- nrow(X)
  n <- ncol(X)
  
  #generate matrix of ranks (optionally keep l2fc's)
  # R <- apply(X, 2, function(x, p) as.integer(rank(x)), p) #original ssgsea
  #balanced ssgsea, with NAs now kept as NA to preserve matrix shape; they're dropped later by order
  if(useranks) R <- apply(X, 2, function(x, p) (as.integer(rank(x, na.last = "keep")) - (sum(!is.na(x))/2 + .5))) else R = X 
  # if(useranks) R <- apply(X, 2, function(x, p) (as.integer(rank(x)) - (p/2 + .5)),p) else R = X #NA's handled incorrectly
  
  #walkplot example
  # if(n == 1 & length(geneSets) == 1){
  #   plotRanking <- order(R[, 1], decreasing = TRUE)
  #   output = walkplot(geneSets[[1]], plotRanking, 1, R, alpha)
  #   return(output)
  # }
  
  #does the summed KS staistic for each sample
  es <- sapply(1:n, function(j, R, geneSets, alpha) {
    geneRanking <- order(R[, j], decreasing = TRUE, na.last = NA) #now removes NA's to ignore them in rndwalk
    es_sample <- NA
    es_sample <- sapply(geneSets, myrndwalk, geneRanking, j, R, alpha)
    unlist(es_sample)
  }, R, geneSets, alpha)
  
  if (length(geneSets) == 1) es <- matrix(es, nrow = 1)
  rownames(es) <- names(geneSets)
  colnames(es) <- colnames(X)
  
  return(es)
}

#GSVA:::rndwalk but able to handle sum(Ralpha) = 0
myrndwalk = function(gSetIdx, geneRanking, j, R, alpha) {
  indicatorFunInsideGeneSet <- match(geneRanking, gSetIdx)
  indicatorFunInsideGeneSet[!is.na(indicatorFunInsideGeneSet)] <- 1
  indicatorFunInsideGeneSet[is.na(indicatorFunInsideGeneSet)] <- 0
  Ralpha = (abs(R[geneRanking, j]) * indicatorFunInsideGeneSet)^alpha
  if(sum(Ralpha) ==0 || sum(!indicatorFunInsideGeneSet) == 0) return(0) #special cases where all or none of gene are in the bag or only one gene with middle rank
  stepCDFinGeneSet <- cumsum(Ralpha)/sum(Ralpha)
  stepCDFoutGeneSet <- cumsum(!indicatorFunInsideGeneSet)/sum(!indicatorFunInsideGeneSet)
  walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet
  return(sum(walkStat))
}

# walkplot= function (gSetIdx, geneRanking, j, R, alpha) 
# {
#   indicatorFunInsideGeneSet <- match(geneRanking, gSetIdx)
#   indicatorFunInsideGeneSet[!is.na(indicatorFunInsideGeneSet)] <- 1
#   indicatorFunInsideGeneSet[is.na(indicatorFunInsideGeneSet)] <- 0
#   stepCDFinGeneSet <- cumsum((abs(R[geneRanking, j]) * indicatorFunInsideGeneSet)^alpha)/sum((abs(R[geneRanking,
#                                                                                                     j]) * indicatorFunInsideGeneSet)^alpha)
#   stepCDFoutGeneSet <- cumsum(!indicatorFunInsideGeneSet)/sum(!indicatorFunInsideGeneSet)
#   walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet
#   return(walkStat)
# }