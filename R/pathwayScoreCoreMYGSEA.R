#' Pathway Score Core - MYGSEA
#' 
#' Computes pathway scores for mygsea.
#' 
#' This function is a parallelized wrapper for MYGSEA, which does the actual
#' scoring. mygsea method uses ranks and normalization, while mygsea_norank
#' method does not use ranks or normalization. Normalization divides final 
#' scores by difference between max and min score. Without normalization,
#' scores from individual samples have no impact on each other. Final
#' pathscoremat is written to disk.
#'
#' @param sk.list Sample keys to use; should correspond to fcmat rownames.
#' @param method Method name to use in file output. "mygsea" or "mygsea_norank"
#' @param pathset Name of pathway set.
#' @param dataset Name of data set.
#' @param fcmat Sample by gene matrix of log2(fold change)'s. Rownames are
#'   sample keys and colnames are genes.
#' @param chem_dict Dataframe with one row per sample key and seven columns:
#'   sample_key, sample_id, conc, time, casrn, name, dsstox_substance_id.
#' @param pathway_data Named ist of gene name vectors. Each element is one 
#'   pathway, defined by the genes it contains.
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
pathwayScoreCoreMYGSEA <- function(sk.list,
                                   method = "mygsea",
                                   pathset="bhrr",
                                   dataset,
                                   fcmat,
                                   chem_dict,
                                   pathway_data,
                                   mc.cores=1,
                                   normalization = T,
                                   useranks = T) {
  
  #change NA to zero, although this is no longer necessary; MYGSEA can handle NA now
  data <- fcmat[sk.list,]
  data[is.na(data)] <- 0

  success <- F
  cat("   start MYGSEA:",dim(data),"\n")
  tryCatch({
    #call MYGSEA for scoring
    if(mc.cores <= 1){
      res <- MYGSEA(X=t(data),
                geneSets =pathway_data,
                min.sz=1,
                max.sz=Inf,
                verbose=T,
                useranks = useranks)
    } else {
      starts = (0:(mc.cores-1)*nrow(data))%/%mc.cores+1
      ends = (1:mc.cores*nrow(data))%/%mc.cores
      datalist = lapply(1:length(starts), function(i){t(data)[,seq(starts[i], ends[i]), drop = F ] } )
      # reslist = mclapply(datalist, MYGSEA, geneSets =pathway_data, min.sz=10, max.sz=Inf, verbose=T, useranks = useranks, 
      #                    mc.cores= mc.cores)

      
      cl = makePSOCKcluster(mc.cores)
      clusterExport(cl, c("myrndwalk"))
      output = clusterEvalQ(cl, library(GSVA))
      reslist = parLapply(cl = cl, X=datalist,fun=MYGSEA,geneSets =pathway_data, min.sz=10, max.sz=Inf, verbose=T, useranks = useranks)
      stopCluster(cl)
      
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
  names(res2) <- c("pathway","sample_key","pathway_score")
  
  #do actual pathway sizes here
  notnamat = !is.na(fcmat)
  trimpdata = lapply(pathway_data, function(x){x[x %in% colnames(fcmat)]})
  sizemat = sapply(trimpdata, function(x){rowSums(notnamat[,x, drop = F])})
  sizemelt = reshape2::melt(sizemat, varnames = c("sample_key", "pathway"), value.name = "size")
  res2$size = sizemelt$size[match(paste0(res2$sample_key, res2$pathway),
                             paste0(sizemelt$sample_key, sizemelt$pathway))]
  
  #construct empty data frame
  sk.list <- res2[,"sample_key"]
  n <- length(sk.list)
  name.list <- c("sample_id","dsstox_substance_id","casrn","name","time","conc","pathset")
  nmat <- as.data.frame(matrix(nrow=n,ncol=length(name.list)))
  colnames(nmat) <- name.list
  cat("   start building output\n")
  
  #add sample identifiers
  dictmatch = match(sk.list,chem_dict$sample_key)
  nmat[,"pathset"] <- pathset
  nmat[,"sample_id"] <- chem_dict$sample_id[dictmatch]
  nmat[,"conc"] <- chem_dict$conc[dictmatch]
  nmat[,"time"] <- chem_dict$time[dictmatch]
  nmat[,"dsstox_substance_id"] <- chem_dict$dsstox_substance_id[dictmatch]
  nmat[,"casrn"] <- chem_dict$casrn[dictmatch]
  nmat[,"name"] <- chem_dict$name[dictmatch]

  #bind scores to ids and trim any unnecessary columns
  pathscoremat <- cbind(nmat,res2)
  cat("   finish building output\n")
  name.list <- c("sample_id","dsstox_substance_id","casrn","name","time","conc","pathset","pathway","size","pathway_score")
  pathscoremat <- pathscoremat[,name.list]
  path.list <- sort(unique(pathscoremat[,"pathway"]))
  
  #write to disk
  file <- paste("output/pathway_score_summary/PATHSCOREMAT_",pathset,"_",dataset,"_",method,".RData",sep="")
  save(pathscoremat,file=file)
  cat("   output written\n")
}
