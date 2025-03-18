#' Signature Score Core - GSEA
#'
#' Computes signature scores for gsea.
#'
#' This function is a parallelized wrapper for gsea, which does the actual
#' scoring. gsea method uses ranks and normalizationby default. Normalization divides final
#' scores by difference between max and min score. Without normalization,
#' scores from individual samples have no impact on each other.
#'
#' @param sk.list Sample keys to use; should correspond to fcmat rownames.
#' @param normfactor = proceed on first 1/normfactor of GSEA data
#' @param sigset Name of signature set.
#' @param fcmat Expects FCMAT2. Sample by gene matrix of log2(fold change)'s. Rownames are
#'   sample keys and colnames are genes.
#' @param chem_dict Expects CHEM_DICT object. Dataframe with one row per sample key and seven columns:
#'   sample_key, sample_id, conc, time, casrn, name, dtxsid.
#' @param signature_data Named list of gene name vectors. Each element is one
#'   signature, defined by the genes it contains.
#' @param mc.cores Number of cores to use.
#' @param normalization normalization = T normalizes final scores.
#' @param useranks useranks = T uses score ranks for weighting; otherwise,
#'   fold changes are used for weights.
#'
#' @import GSVA
#' @importFrom parallel makePSOCKcluster clusterExport clusterEvalQ parLapply stopCluster
#' @importFrom reshape2 melt
#'
#' @return Dataframe with one row per chemical/conc/signature combination. Columns
#'   are: sample_id, dtxsid, casrn, name, time, conc, sigset, signature, size, signature_score
#' @export signatureScoreCoreGSEA
signatureScoreCoreGSEA <- function(sk.list,
                                   normfactor=7500,
                                   sigset,
                                   fcmat,
                                   chem_dict,
                                   signature_data,
                                   mc.cores=1,
                                   normalization = T,
                                   useranks = T) {

  printCurrentFunction(paste(sigset))
  #change NA to zero, although this is no longer necessary; gsea can handle NA now
  data <- fcmat[sk.list,]
  data[is.na(data)] <- 0

  success <- F
  cat("   start GSEA:",dim(data),":",mc.cores,"\n")
  tryCatch({
    #call GSEA for scoring
    if(mc.cores <= 1){
      res <- GSEA(X=t(data),
                  geneSets =signature_data,
                  min.sz=1,
                  max.sz=Inf,
                  verbose=F,
                  useranks = useranks)
    } else {
      starts = (0:(mc.cores-1)*nrow(data))%/%mc.cores+1
      ends = (1:mc.cores*nrow(data))%/%mc.cores
      cat("   make datalist\n")
      datalist = lapply(1:length(starts), function(i){t(data)[,seq(starts[i], ends[i]), drop = F ] } )
      # reslist = mclapply(datalist, GSEA, geneSets =signature_data, min.sz=1, max.sz=Inf, verbose=T, useranks = useranks,
      #                    mc.cores= mc.cores)


      cat("   make cluster\n")
      cl = makePSOCKcluster(mc.cores)
      cat("   export cluster\n")
      clusterExport(cl, c("myrndwalk"))
      cat("   eval cluster\n")
      output = clusterEvalQ(cl, library(GSVA))
      cat("   make reslist\n")
      # [LJE 5/13/24]
      # Changed min.sz=1 (was 10) in the line below to match how GSEA function is invoked when calling with mc.cores = 1
      # This was enforcing minimum signature size of at least 10 even when minsigsz was set to a lower number in the driver function
      # Based on the documentation in GSEA.R that min.sz is a deprecated parameter, it probably should not be set here at all
      # Filtering by signature size is now handled completely before this function is called.
      reslist = parLapply(cl = cl, X=datalist,fun=GSEA,geneSets =signature_data, min.sz=1, max.sz=Inf, verbose=F, useranks = useranks)
      cat("   stop cluster\n")
      stopCluster(cl)
      cat("   reduce\n")
      res = Reduce(cbind, reslist)
    }
    if (normalization) {
      #nfactor = (range(res)[2] - range(res)[1])
      # this old value made values for a specific signature dependent on the other signatures
      nfactor = normfactor
      res = res/nfactor
    }
    success <- T
  }, warning = function(w) {
    print(w)
    success <- T
  }, error = function(e) {
    print(e)
    success <- F
  })
  cat("   finish GSEA: ",dim(res),"\n")

  #melt output
  res2 <- melt(res,as.is=T)
  names(res2) <- c("signature","sample_key","signature_score")

  #do actual signature sizes here
  notnamat = !is.na(fcmat)
  trimpdata = lapply(signature_data, function(x){x[x %in% colnames(fcmat)]})
  sizemat = sapply(trimpdata, function(x){rowSums(notnamat[,x, drop = F])})
  sizemelt = melt(sizemat, varnames = c("sample_key", "signature"), value.name = "size")
  res2$size = sizemelt$size[match(paste0(res2$sample_key, res2$signature),
                                  paste0(sizemelt$sample_key, sizemelt$signature))]

  #construct empty data frame
  sk.list <- res2[,"sample_key"]
  n <- length(sk.list)
  name.list <- c("sample_id","dtxsid","casrn","name","time","conc","sigset")
  nmat <- as.data.frame(matrix(nrow=n,ncol=length(name.list)))
  colnames(nmat) <- name.list
  cat("   start building output\n")

  #add sample identifiers
  dictmatch = match(sk.list,chem_dict$sample_key)
  nmat[,"sigset"] <- sigset
  nmat[,"sample_id"] <- chem_dict$sample_id[dictmatch]
  nmat[,"conc"] <- chem_dict$conc[dictmatch]
  nmat[,"time"] <- chem_dict$time[dictmatch]
  nmat[,"dtxsid"] <- chem_dict$dtxsid[dictmatch]
  nmat[,"casrn"] <- chem_dict$casrn[dictmatch]
  nmat[,"name"] <- chem_dict$name[dictmatch]

  #bind scores to ids and trim any unnecessary columns
  signaturescoremat <- cbind(nmat,res2)
  cat("   finish building output\n")
  if(is.element("time",names(signaturescoremat)))
    name.list <- c("sample_id","dtxsid","casrn","name","time","conc","sigset","signature","size","signature_score")
  else
    name.list <- c("sample_id","dtxsid","casrn","name","conc","sigset","signature","size","signature_score")

  signaturescoremat <- signaturescoremat[,name.list]
  signature.list <- sort(unique(signaturescoremat[,"signature"]))

  return(signaturescoremat)


}
