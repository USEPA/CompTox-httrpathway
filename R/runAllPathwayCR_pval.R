#' Run All Pathway Concentration Response (P-Value)
#' 
#' Driver for pathway scoring and concentration response (CR).
#' 
#' CR requires pathway scores to have already been computed for a nullset.
#' randomdata() can generate a nullset, and this function can compute pathway
#' scores for it by setting dataset = nullset and nullset = NULL. Pathway
#' scores are written to disk in output/pathway_score_summary/. CR results
#' are written to disk in output/pathway/conc_resp_summary/.
#'
#' @param basedir Folder that stores FCMAT2 and CHEM_DICT files.
#' @param dataset Name of data set.
#' @param pathset Name of pathway set.
#' @param method Name of pathway scoring method.
#' @param minpathsize Minimum pathway size. 
#' @param conthits conthits = T uses continous hitcall; conthits = F uses discrete
#'   hitcalls.
#' @param nullset Name of null dataset. Set nullset = NULL to skip CR.
#' @param do.plot do.plot = T generates a CR plot for every sample/pathway
#'   combination.
#' @param pval P-value to use for noise estimation.
#' @param mc.cores Vector with two values: number of cores to use for pathway
#'   scoring and number of cores to use for CR. CR can usually handle the maximum
#'   number, but gsva scoring might require a smaller number to avoid memory
#'   overflow.
#' @param fitmodels Vector of model names to run conc/resp with. "cnst" should 
#'   always be chosen.
#'   
#' @return No output.
#' @export
runAllPathwayCR_pval = function(basedir="input/fcdata/",dataset="arer",pathset="bhrr", 
                                method = "fc", minpathsize = 10, conthits = T,
                                nullset = "arer_RAND125",do.plot = T, pval = .05, 
                                mc.cores = c(39,39), 
                                fitmodels = c("cnst", "hill", "gnls", "poly1", "poly2", 
                                              "pow", "exp2", "exp3", "exp4", "exp5")){
  
  # load fcmat, chemdict   
  file <- paste0(basedir,"FCMAT2_",dataset,".RData")
  load(file)
  file <- paste0(basedir,"CHEM_DICT_",dataset,".RData",sep="")
  load(file)
  rownames(CHEM_DICT) <- CHEM_DICT[,"sample_key"]
  
  #run pathway scores
  pathwayScore(FCMAT2, CHEM_DICT, pathset,dataset,method=method,mc.cores=mc.cores[1], minpathsize = minpathsize)
  if(is.null(nullset)) return()
  
  #run conc/resp
  pathwayConcResp_pval(pathset,dataset,method=method, nullset = nullset, mc.cores=mc.cores[2], do.plot = do.plot,
                       to.file=T, pval = pval, minpathsize = minpathsize, conthits = conthits,
                       fitmodels = fitmodels)

}

