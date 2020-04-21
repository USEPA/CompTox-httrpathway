#' Signature Score Core - MYGSEA
#'
#' Computes signature scores for mygsea.
#'
#' This function is a parallelized wrapper for MYGSEA, which does the actual
#' scoring. mygsea method uses ranks and normalization, while mygsea_norank
#' method does not use ranks or normalization. Normalization divides final
#' scores by difference between max and min score. Without normalization,
#' scores from individual samples have no impact on each other. Final
#' signaturescoremat is written to disk.
#'
#' @param sk.list Sample keys to use; should correspond to fcmat rownames.
#' @param method Method name to use in file output. "mygsea" or "mygsea_norank"
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
#' @param normalization normalization = T normalizes final scores.
#' @param useranks useranks = T uses score ranks for weighting; otherwise,
#'   fold changes are used for weights.
#'
#' @import openxlsx
#' @import GSVA
#' @import parallel
#'
#' @return No output.
#' @export
signatureScoreCoreMYGSEA <- function(sk.list,
                                     method = "mygsea",
                                     sigset,
                                     dataset,
                                     fcmat,
                                     chem_dict,
                                     signature_data,
                                     mc.cores=1,
                                     normalization = T,
                                     useranks = T) {

  printCurrentFunction(paste(dataset,sigset,method))
  #change NA to zero, although this is no longer necessary; MYGSEA can handle NA now
  data <- fcmat[sk.list,]
  data[is.na(data)] <- 0

  success <- F
  cat("   start MYGSEA:",dim(data),":",mc.cores,"\n")
  tryCatch({
    #call MYGSEA for scoring
    if(mc.cores <= 1){
      res <- MYGSEA(X=t(data),
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
      # reslist = mclapply(datalist, MYGSEA, geneSets =signature_data, min.sz=10, max.sz=Inf, verbose=T, useranks = useranks,
      #                    mc.cores= mc.cores)


      cat("   make cluster\n")
      cl = makePSOCKcluster(mc.cores)
      cat("   export cluster\n")
      clusterExport(cl, c("myrndwalk"))
      cat("   eval cluster\n")
      output = clusterEvalQ(cl, library(GSVA))
      cat("   make reslist\n")
      reslist = parLapply(cl = cl, X=datalist,fun=MYGSEA,geneSets =signature_data, min.sz=10, max.sz=Inf, verbose=F, useranks = useranks)
      cat("   stop cluster\n")
      stopCluster(cl)
      cat("   reduce\n")
      res = Reduce(cbind, reslist)
    }
    if (normalization) {
      nfactor = (range(res)[2] - range(res)[1])
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
  cat("   finish mygsea: ",dim(res),"\n")

  #melt output
  res2 <- reshape2::melt(res,as.is=T)
  names(res2) <- c("signature","sample_key","signature_score")

  #do actual signature sizes here
  notnamat = !is.na(fcmat)
  trimpdata = lapply(signature_data, function(x){x[x %in% colnames(fcmat)]})
  sizemat = sapply(trimpdata, function(x){rowSums(notnamat[,x, drop = F])})
  sizemelt = reshape2::melt(sizemat, varnames = c("sample_key", "signature"), value.name = "size")
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

  #write to disk
  file <- paste("../output/signature_score_summary/signaturescoremat_",sigset,"_",dataset,"_",method,"_directional.RData",sep="")
  cat("   savging file ... ",file,"\n")
  save(signaturescoremat,file=file)
  cat("   output written\n")
}
