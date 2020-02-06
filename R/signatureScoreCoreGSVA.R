#' Signature Score Core - GSVA
#'
#' Computes GSVA signature scores.
#'
#' This function is a wrapper for GSVA with Gaussian cdf kernels. signaturescoremat
#' output is saved directly to disk.
#'
#'
#' @param sk.list Sample keys to use; should correspond to fcmat rownames.
#' @param sigset Name of signature set.
#' @param dataset Name of data set.
#' @param fcmat Sample by gene matrix of log2(fold change)'s. Rownames are
#'   sample keys and colnames are genes.
#' @param chem_dict Dataframe with one row per sample key and seven columns:
#'   sample_key, sample_id, conc, time, casrn, name, dtxsid.
#' @param signature_data Named ist of gene name vectors. Each element is one
#'   signature, defined by the genes it contains.
#' @param mc.cores Number of cores to use. Parallelization is performed
#'   by gsva itself.
#'
#' @import openxlsx
#' @import GSVA
#' @import parallel
#'
#' @return No output.
#' @export
signatureScoreCoreGSVA <- function(sk.list,
                                   sigset="FILTERED",
                                   dataset,
                                   fcmat,
                                   chem_dict,
                                   signature_data,
                                   mc.cores=1) {

  #keep only rows in sk.list and change NA's to 0's because gsva can't handle NA
  #This should cause NA's to cluster in the middle of the K-S statistic and have minimal effect
  #on the score. Ideally, GSVA would handle NA's internally.
  data <- fcmat[sk.list,]
  data[is.na(data)] <- 0
  data = data[,!apply(data, 2, function(x){all(x == x[1])})] #remove constant columns

  success <- F
  cat("   start gsva:",dim(data),"\n")
  #call gsva: converts gene by sample matrix to signature by sample matrix
  tryCatch({
    res <- gsva(expr=t(data),
                gset.idx.list=signature_data,
                method="gsva",
                min.sz=1,
                max.sz=Inf,
                parallel.sz=mc.cores,
                verbose=T,
                kcdf = "Gaussian"
                )
    success <- T
  }, warning = function(w) {
    print(w)
    success <- T
  }, error = function(e) {
    print(e)
    success <- F
  })
  cat("   finish gsva: ",dim(res),"\n")

  #melt the gsva output
  res2 <- reshape2::melt(res,as.is=T)
  names(res2) <- c("signature","sample_key","signature_score")

  #do actual signature sizes here
  notnamat = !is.na(fcmat) #get nonmissings
  trimpdata = lapply(signature_data, function(x){x[x %in% colnames(fcmat)]}) #remove genes not in fcmat from signature_data
  sizemat = sapply(trimpdata, function(x){rowSums(notnamat[,x, drop = F])}) #count nonmissings and out skey x signature matrix
  sizemelt = reshape2::melt(sizemat, varnames = c("sample_key", "signature"), value.name = "size") #melt it so we can do a match
  res2$size = sizemelt$size[match(paste0(res2$sample_key, res2$signature),
                                  paste0(sizemelt$sample_key, sizemelt$signature))] #match pasted samplekey and signature

  #create skeleton of output data frame
  sk.list <- res2[,"sample_key"]
  n <- length(sk.list)
  name.list <- c("sample_id","dtxsid","casrn","name","time","conc","sigset")
  nmat <- as.data.frame(matrix(nrow=n,ncol=length(name.list)))
  colnames(nmat) <- name.list
  cat("   start building output\n")

  #fill in identifiers using chem_dict
  dictmatch = match(sk.list,chem_dict$sample_key)
  nmat[,"sigset"] <- sigset
  nmat[,"sample_id"] <- chem_dict$sample_id[dictmatch]
  nmat[,"conc"] <- chem_dict$conc[dictmatch]
  nmat[,"time"] <- chem_dict$time[dictmatch]
  nmat[,"dtxsid"] <- chem_dict$dtxsid[dictmatch]
  nmat[,"casrn"] <- chem_dict$casrn[dictmatch]
  nmat[,"name"] <- chem_dict$name[dictmatch]

  #attach melted scores and strip unusued columns
  signaturescoremat <- cbind(nmat,res2)
  cat("   finish building output\n")
  name.list <- c("sample_id","dtxsid","casrn","name","time","conc","sigset","signature","size","signature_score")
  signaturescoremat <- signaturescoremat[,name.list]
  signature.list <- sort(unique(signaturescoremat[,"signature"]))

  #write output to disk
  method <- "gsva"
  file <- paste("../output/signature_score_summary/signaturescoremat_",sigset,"_",dataset,"_",method,"_directional.RData",sep="")
  save(signaturescoremat,file=file)
  cat("   output written\n")
}
