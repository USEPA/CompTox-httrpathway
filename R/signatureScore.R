#' Signature Score
#'
#' Computes and saves signature scores.
#'
#' signatureScore is a driver for various scoring methods. The three that are
#' currently available are "gsva", "gsea", "fc", and "gsea_norank" (a version
#' of gsea that uses fold changes instead of ranks as weights). Deprecated
#' methods include the Fisher method and gsvae (gsva with empirical cdfs).
#' Beware running out of memory on large runs with gsva, Linux, and many cores.
#' Signature size is counted according to number of genes in the signature that are
#' also in the column names of FCMAT2. However, each method performs a more
#' rigorous size count internally that accounts for missing values and adds this
#' to the output. This minsigsize is enforced when running signatureConcResp_pval.
#'
#' @param FCMAT2 Sample by gene matrix of log2(fold change)'s. Rownames are
#'   sample keys and colnames are genes.
#' @param CHEM_DICT Dataframe with one row per sample key and seven columns:
#'   sample_key, sample_id, conc, time, casrn, name, dtxsid.
#' @param sigset Name of signature set.
#' @param sigcatalog Name of the signature catalog file
#' @param dataset Name of data set.
#' @param method Signature scoring method in c("fc", "gsva", "gsea")
#' @param normfactor Value passed ot the plotting code to scale the y values
#' @param mc.cores Number of cores to use.
#' @param minsigsize Minimum allowed signature size BEFORE accounting for
#'   missing values.
#'
#' @import data.table
#' @import parallel
#' @import openxlsx
#'
#' @return No output.
#' @export
signatureScore <- function(FCMAT2,
                           CHEM_DICT,
                           sigset,
                           sigcatalog,
                           dataset,
                           method,
                           normfactor=7500,
                           mc.cores=1,
                           minsigsize = 10) {

  printCurrentFunction(paste(dataset,sigset,method))
  starttime = proc.time()
  cat("  signatureScore: create directory\n")
  dir.create("../output/signature_score_summary/", showWarnings = F)

  cat("  signatureScore: Start nonempties\n")

  #get rid of columns filled with missing values
  nonempties = apply(FCMAT2,2,function(x){sum(!is.na(x))})
  FCMAT2 = FCMAT2[,nonempties > 0]

  cat("  signatureScore: load signature data\n")
  if(sigset!="wgcna") {
    file = paste0("../input/signatures/signatureDB_genelists.RData")
    print(file)
    load(file=file)
  }
  else {
    name <- str_replace(sigcatalog,"_catalog","")
    file = paste0("../input/signatures/",name,".RData")
    print(file)
    load(file=file)
    genelists <- sigdb
  }

  #file <- paste0("../input/signatures/signatureDB_genelists.RData")
  #cat("   ",file,"\n")
  #load(file) #genelists
  cat("  signatureScore: signatureCatalogLoader\n")
  catalog <- signatureCatalogLoader(sigset,sigcatalog)

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
    file <- paste0("../output/deleted signatures ",dataset," ",sigset,".xlsx")
    write.xlsx(slost,file)
    browser()
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
                             sigset = sigset,
                             dataset = dataset,
                             chem_dict = CHEM_DICT,
                             signature_data = signature_data )
      stopCluster(cl)
      cat("  signatureScore: finish signatureScoreCoreFC parallel\n")
      #reform output
      signaturescoremat = as.data.frame(rbindlist(pscorelist))
    } else {
      cat("  signatureScore: Start signatureScoreCoreFC single\n")
      signaturescoremat = signatureScoreCoreFC(fcdata = FCMAT2, sigset = sigset,dataset = dataset,
                                          chem_dict = CHEM_DICT,signature_data = signature_data )
      cat("  signatureScore: Finish signatureScoreCoreFC single\n")
    }

    cat("  signatureScore: save signaturescoremat\n")
    file <- paste0("../output/signature_score_summary/signaturescoremat_",sigset,"_",dataset,"_",method,"_directional.RData")
    cat("   ",file,"\n")
    #browser()
    save(signaturescoremat,file=file)
  }
  #call gsva scoring
  if(method=="gsva") {
    signatureScoreCoreGSVA(sk.list,
                           sigset=sigset,
                           dataset=dataset,
                           fcmat=FCMAT2,
                           chem_dict=CHEM_DICT,
                           signature_data=signature_data,
                           mc.cores=mc.cores)
  }
  #call gsea scoring
  if(method=="gsea") {
    signatureScoreCoregsea(sk.list,
                             method = "gsea",
                             normfactor=normfactor,
                             sigset=sigset,
                             dataset=dataset,
                             fcmat=FCMAT2,
                             chem_dict=CHEM_DICT,
                             signature_data=signature_data,
                             mc.cores=mc.cores,
                             normalization = T,
                             useranks = T)
  }
  #call gsea scoring without ranks or normalization
  if(method=="gsea_norank") {
    signatureScoreCoregsea(sk.list,
                             method = "gsea_norank",
                             sigset=sigset,
                             dataset=dataset,
                             fcmat=FCMAT2,
                             chem_dict=CHEM_DICT,
                             signature_data=signature_data,
                             mc.cores=mc.cores,
                             normalization = F,
                             useranks = F)
  }

  cat("\n  signatureScore: Time taken:",proc.time() - starttime, "\n")

}
