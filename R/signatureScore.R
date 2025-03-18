#' Signature Score
#'
#' Computes and saves signature scores.
#'
#' signatureScore is a driver for various scoring methods. The three that are
#' currently available are "gsva", "gsea", "fc".
#' Beware running out of memory on large runs with gsva, Linux, and many cores -- ensure your system has enough memory allocated depending on data size.
#' Signature size is counted according to number of genes in the signature that are
#' also in the column names of FCMAT2. However, each method performs a more
#' rigorous size count internally that accounts for missing values and adds this
#' to the output.
#'
#' @param FCMAT2 Sample by gene matrix of log2(fold change)'s. Rownames are
#'   sample keys and colnames are genes.
#' @param CHEM_DICT Dataframe with one row per sample key and seven columns:
#'   sample_key, sample_id, conc, time, casrn, name, dtxsid.
#' @param sigset Name of signature set.
#' @param sigcatalog full path to signature catalog xlsx file; default is repo version
#' @param method Signature scoring method in c("fc", "gsva", "gsea")
#' @param normfactor Value passed to the plotting code to scale the y values
#' @param mc.cores Number of cores to use.
#' @param minsigsize Minimum allowed signature size BEFORE accounting for
#'   missing values.
#' @param sigdbgenelist full path to signature DB gene list file; default is repo version
#'
#' @import data.table
#' @importFrom parallel makePSOCKcluster clusterExport parLapply stopCluster
#' @importFrom openxlsx write.xlsx
#' @importFrom stringr str_replace
#'
#' @return Returns data frame of signature scores
#' @export signatureScore
signatureScore <- function(FCMAT2,
                           CHEM_DICT,
                           sigset,
                           sigcatalog = "../inst/extdata/signatureDB_master_catalog_2022-05-16.xlsx",
                           method,
                           normfactor=7500,
                           mc.cores=1,
                           minsigsize = 10,
                           sigdbgenelist = "../inst/extdata/signatureDB_genelists.RDS") {

  printCurrentFunction(paste(sigset,method))
  starttime = proc.time()

  cat("  signatureScore: Start nonempties\n")

  #get rid of columns filled with missing values
  nonempties = apply(FCMAT2,2,function(x){sum(!is.na(x))})
  FCMAT2 = FCMAT2[,nonempties > 0]

  cat("  signatureScore: load signature data\n")
  print(sigdbgenelist)
  genelists <- readRDS(sigdbgenelist)


  cat("  signatureScore: signatureCatalogLoader\n")
  catalog <- signatureCatalogLoader(sigset,sigcatalog,sigdbgenelist)

  catalog <- catalog[is.element(catalog$signature,names(genelists)),]
  signature_data <- genelists[catalog$signature]

  #sk.list could be used to choose a data subset, but here does nothing
  sk.list <-as.matrix(rownames(FCMAT2))

  #Enforce minimum signature size
  plengths = sapply(signature_data, function(x){sum(x %in% colnames(FCMAT2))})
  cat("  Removing", sum(plengths < minsigsize), "signatures under min size", minsigsize, ".", sum(plengths >= minsigsize),
      "signatures remaining.\n")
  signature_data0 <- signature_data
  signature_data = signature_data[plengths >= minsigsize]

  sig.list.0 <- sort(unique(names(signature_data0)))
  sig.list <- sort(unique(names(signature_data)))
  slost <- sig.list.0[!is.element(sig.list.0,sig.list)]
  if(length(slost)>0) {
    for(i in 1:length(slost)) cat("   ",slost[i],"\n")
  }
  cat("  signatureScore: now run the inner functions depending on method\n")
  if(method=="fc") {
    if(mc.cores > 1){
      #split fcmat into mc.cores matrices and run them in parallel
      splitseq = sort(rep(1:mc.cores, length.out = nrow(FCMAT2)))
      fclist = split.data.frame(FCMAT2,splitseq, drop = F)

      cl = makePSOCKcluster(mc.cores)
      clusterExport(cl, c("signatureScoreCoreFC"))
      cat("  signatureScore: start signatureScoreCoreFC parallel\n")
      pscorelist = parLapply(cl = cl, X=fclist, fun=signatureScoreCoreFC,
                             fcmat = FCMAT2,
                             sigset = sigset,
                             chem_dict = CHEM_DICT,
                             signature_data = signature_data )
      stopCluster(cl)
      cat("  signatureScore: finish signatureScoreCoreFC parallel\n")
      #reform output
      score = as.data.frame(rbindlist(pscorelist))
    } else {
      cat("  signatureScore: Start signatureScoreCoreFC single\n")
      score <- signatureScoreCoreFC(fcmat = FCMAT2, sigset = sigset,
                                    chem_dict = CHEM_DICT,signature_data = signature_data )
      cat("  signatureScore: Finish signatureScoreCoreFC single\n")
    }
  }
  #call gsva scoring
  if(method=="gsva") {
    score <- signatureScoreCoreGSVA(sk.list,
                                    sigset=sigset,
                                    fcmat=FCMAT2,
                                    chem_dict=CHEM_DICT,
                                    signature_data=signature_data,
                                    mc.cores=mc.cores)
  }
  #call gsea scoring
  if(method=="gsea") {
    score <- signatureScoreCoreGSEA(sk.list,
                                    normfactor=normfactor,
                                    sigset=sigset,
                                    fcmat=FCMAT2,
                                    chem_dict=CHEM_DICT,
                                    signature_data=signature_data,
                                    mc.cores=mc.cores,
                                    normalization = T,
                                    useranks = T)
  }
  #call gsea scoring without ranks or normalization
  if(method=="gsea_norank") {
    score <- signatureScoreCoreGSEA(sk.list,
                                    sigset=sigset,
                                    fcmat=FCMAT2,
                                    chem_dict=CHEM_DICT,
                                    signature_data=signature_data,
                                    mc.cores=mc.cores,
                                    normalization = F,
                                    useranks = F)
  }

  cat("\n  signatureScore: Time taken:",proc.time() - starttime, "\n")
  return(score)

}
