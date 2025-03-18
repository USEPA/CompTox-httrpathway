#' Randomized Null Data
#'
#' Generate randomized null data based on actual data.
#'
#' Randomization is performed by sampling the quantile function for each gene in the actual data.
#' The null set will have roughly the same distribution of values for each gene in the actual data.
#'
#' @param nchem Number of null chemicals. Number of null samples is approximately
#'   eight times this value.
#' @param seed Random seed.
#' @param maxconc Only use concentrations less than maxconc, default 1000000
#' @param nlowconc If not NULL, only include the lowest nlowconc concentrations for each chemical
#' @param dtxsid.exclude dtxsids to exclude, default NULL
#' @param FCMAT2 Sample by gene matrix of log2(fold change)'s. Rownames are sample keys and colnames are genes.
#' @param CHEM_DICT Data frame with one row per sample key and seven columns:
#'   sample_key, sample_id, conc, time, casrn, name, dtxsid.
#' @param output_dir full path to directory/file where the CHEM_DICT/FCMAT2 data frames should be saved
#' @param writing_to_disk boolean, TRUE if the resulting CHEM_DICT/FCMAT2 data frames are to be saved to disk
#' @importFrom stats quantile
#' @importFrom stats runif
#' @return list of the randomized FCMAT2 and CHEM_DICT data frame
#' @export randomdata
randomdata <- function(nchem = 1000,
                      seed = 12345,
                      maxconc=1000000,
                      nlowconc=2,
                      dtxsid.exclude=NULL,
                      FCMAT2,
                      CHEM_DICT,
                      output_dir="../data/",
                      writing_to_disk=FALSE){
  printCurrentFunction()
  set.seed(seed)

  chem_dict_filter = CHEM_DICT
  if(!is.null(nlowconc)) {
    filter = vector(length=nrow(chem_dict_filter),mode="integer")
    filter[] = 0
    slist = unique(chem_dict_filter$sample_id)
    cdf2 = NULL
    for(sid in slist) {
      temp = chem_dict_filter[is.element(chem_dict_filter$sample_id,sid),]
      clist =sort(unique(temp$conc))
      temp = temp[temp$conc<=clist[nlowconc],]
      cdf2 = rbind(cdf2,temp)
    }
    chem_dict_filter = cdf2
  }
  chem_dict_filter = chem_dict_filter[chem_dict_filter$conc<maxconc,]

  if(!is.null(dtxsid.exclude)) chem_dict_filter = chem_dict_filter[!is.element(chem_dict_filter$dtxsid,dtxsid.exclude),]
  chem_dict_0 = CHEM_DICT

   #get all the concentration vectors and sample nchem of them with replacement
  concpats = lapply(CHEM_DICT$sample_id, function(x){CHEM_DICT$conc[CHEM_DICT$sample_id == x]}) #get concentration patterns

  #  nconcs = sapply(concpats, function(x) {length(unique(x))})
  nconcs = sapply(concpats, function(x) {length(x)})
  concpats = concpats[nconcs >= 4] #don't use concs with <4 values
  concpats = concpats[nconcs <= 8] #don't use concs with >8 values
  concpats = concpats[sample(1:length(concpats), nchem, replace = T)]
  concs = unlist(concpats)
  clens = sapply(concpats[1:nchem],length)
  n = length(concs)
  sk.list = chem_dict_filter$sample_key
  FCMAT2 = FCMAT2[sk.list,]

  fcmat2_0 <- FCMAT2
  cat("  dim:",dim(FCMAT2),"\n")

  FCMAT2 = apply(FCMAT2, 2, function(x){
    #if all are NA, keep it that way
    if(sum(!is.na(x)) == 0){
      out = rep(NA, n)
    } else if(sum(!is.na(x)) == 1){
      #if only one nonmissing, we will insert at least one (more if n > length) of this value in a random location
      out = rep(x[!is.na(x)], n)
    } else {
      #sample uniform distribution and use quantile function to map onto probability
      #distribution of data
      out = quantile(x, runif(n), na.rm = T)
    }
    #read missing values randomly in same percentage as original
    percmissing = sum(is.na(x))/length(x)
    nanum = floor(percmissing*n)
    out[sample(n, nanum)] = NA

    # #if you want an example
    # plot(density(x[!is.na(x)]), ylab = "Probability Density", xlab = "l2fc", main = "Gene: AAAS")
    # points(density(out[!is.na(out)]), type = "l", col = "red")
    # legend("topleft", legend = c("Empirical", paste0(n, " random points")), col = c("black", "red"), lty = c(1,1))

    return(out)
  })

  #prepare CHEM_DICT and save to disk
  CHEM_DICT = data.frame(sample_key = rep(NA,n), sample_id = rep(paste0("Randomid", 1:nchem),times = clens),
                         conc = concs, time = rep(0,n), casrn = rep(paste0("Randomcasrn", 1:nchem), times = clens),
                         name = rep(paste0("Randomname", 1:nchem), times = clens),
                         dtxsid = rep(paste0("Randomdtx", 1:nchem), times = clens),
                         stringsAsFactors = F)
  CHEM_DICT$sample_key = paste(CHEM_DICT$sample_id, CHEM_DICT$conc, sep = "_")
  if (writing_to_disk==TRUE){
    saveRDS(CHEM_DICT, paste0(output_dir,"CHEM_DICT_RAND",nchem,".RDS"))
  }

  #save fcmat to disk
  rownames(FCMAT2) = CHEM_DICT$sample_key
  if (writing_to_disk==TRUE){
    saveRDS(FCMAT2, paste0(output_dir,"FCMAT2_RAND",nchem,".RDS"))
  }
  return(list(FCMAT2,CHEM_DICT))
}
