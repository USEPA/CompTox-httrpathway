#' Run All Replicate Chemical Concentration Response
#'
#' Runs signature scoring and concentration response for replicate chemicals.
#'
#' This function has hard-coded dataset names for the replicates. For each
#' replicate, it computes signature scores, generates a null dataset, runs
#' signature scores for the null dataset, and then runs concentration-response
#' on the actual data. Pathway scores and CR are written to disk.
#'
#' @param basedir Folder that the FCMAT2's are stored in.
#' @param pathset Name of signature set.
#' @param method Name of signature scoring method.
#' @param minpathsize Minimum signature size.
#' @param do.plot do.plot = T generates plots for every chemical/signature/replicate
#'   combination. Adds a significant amount to the runtime.
#' @param pval P-value to use for noise estimation.
#' @param mc.cores Vector with two values: number of cores to use for signature
#'   scoring and number of cores to use for CR. CR can usually handle the maximum
#'   number, but gsva scoring might require a smaller number to avoid memory
#'   overflow.
#' @param conthits conthits = T uses continuous hitcalls. Continuous hitcalls are
#'   a prerequisitie for using repChemPathwayPlot().
#' @param nchem Number of null chemicals to use. The number of null samples is
#'   approximately eight times this value, so nchem = 125 generates ~1000 null
#'   samples.
#'
#' @return No output.
#' @export
runAllRepChemCR = function(basedir="../input/fcdata/",pathset="bhrr", method = "fc", minpathsize = 10,
                           do.plot = F, pval = .05, mc.cores = c(39,39), conthits = T, nchem = 125){

  #hard-coded dataset names
  studys = c("ph1_", "pilot_")
  floors = c("5","10")
  pes = c("0","1")
  methods = c("none", "normal", "normalold", "apeglm", "ashr")
  combos = expand.grid(studys, floors,pes,methods, stringsAsFactors = F)
  datanames = apply(combos,1,function(x){paste0(c(x,"_gene"),collapse = "")})

  for(dataset in datanames){
    print(dataset)
    #first signature scoring
    file <- paste0(basedir,"FCMAT2_",dataset,".RData")
    load(file)
    file <- paste0(basedir,"CHEM_DICT_",dataset,".RData",sep="")
    load(file)
    rownames(CHEM_DICT) <- CHEM_DICT[,"sample_key"]
    signatureScore(FCMAT2, CHEM_DICT, pathset,dataset,method=method,mc.cores=mc.cores[1], minpathsize = minpathsize)

    #create randomized null data (only need to run this the first time)
    randomdata(dataset=dataset, nchem = nchem, seed = 12345)
    #score randomized data
    nullname = paste0(dataset,"_RAND", nchem)
    file <- paste0(basedir,"FCMAT2_",nullname,".RData")
    load(file)
    file <- paste0(basedir,"CHEM_DICT_",nullname,".RData",sep="")
    load(file)
    rownames(CHEM_DICT) <- CHEM_DICT[,"sample_key"]
    signatureScore(FCMAT2, CHEM_DICT, pathset,nullname,method=method,mc.cores=mc.cores[1], minpathsize = minpathsize)

    #regular CR
    signatureConcResp_pval(pathset,dataset,method=method, nullset = nullname, mc.cores=mc.cores[2], do.plot = do.plot,
                         to.file=T, pval = pval, conthits= conthits, minpathsize = minpathsize)

    # #NULL CR
    # signatureConcResp_pval(pathset,dataset = nullname,method=method, nullset = nullname, mc.cores=mc.cores[2], do.plot = do.plot,
    #                      to.file=T, pval = pval, conthits= conthits, minpathsize = minpathsize)

  }

}
