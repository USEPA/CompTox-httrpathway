#' Pathway Score
#' 
#' Computes and saves pathway scores.
#' 
#' pathwayScore is a driver for various scoring methods. The three that are
#' currently available are "gsva", "mygsea", "fc", and "mygsea_norank" (a version
#' of mygsea that uses fold changes instead of ranks as weights). Deprecated
#' methods include the Fisher method and gsvae (gsva with empirical cdfs). 
#' Beware running out of memory on large runs with gsva, Linux, and many cores.
#' Pathway size is counted according to number of genes in the pathway that are 
#' also in the column names of FCMAT2. However, each method performs a more
#' rigorous size count internally that accounts for missing values and adds this
#' to the output. This minpathsize is enforced when running pathwayConcResp_pval.
#'
#' @param FCMAT2 Sample by gene matrix of log2(fold change)'s. Rownames are
#'   sample keys and colnames are genes.
#' @param CHEM_DICT Dataframe with one row per sample key and seven columns:
#'   sample_key, sample_id, conc, time, casrn, name, dsstox_substance_id.
#' @param pathset Name of pathway set.
#' @param dataset Name of data set.
#' @param method Name of desired scoring method.
#' @param mc.cores Number of cores to use.
#' @param minpathsize Minimum allowed pathway size BEFORE accounting for
#'   missing values. 
#' 
#' @import data.table
#' @import parallel
#' @import openxlsx
#'
#' @return No output.
#' @export
pathwayScore <- function(FCMAT2,
                         CHEM_DICT,
                         pathset="bhrr",
                         dataset="PlateEffect", 
                         method="10chems",
                         mc.cores=1,
                         minpathsize = 10) {
  
  starttime = proc.time()
  dir.create("output/pathway_score_summary/", showWarnings = F)
  
  #get rid of columns filled with missing values
  nonempties = apply(FCMAT2,2,function(x){sum(!is.na(x))})
  FCMAT2 = FCMAT2[,nonempties > 0]
  
  #load pathway data
  file <- paste0("input/processed_pathway_data/PATHWAY_LIST_",pathset,".RData")
  load(file) #pathway_data
  
  #sk.list could be used to choose a data subset, but here does nothing
  sk.list <-as.matrix(rownames(FCMAT2))
  
  #Enforce minimum pathway size 
  plengths = sapply(pathway_data, function(x){sum(x %in% colnames(FCMAT2))})
  cat("Removing", sum(plengths < minpathsize), "pathways under min size", minpathsize, ".", sum(plengths >= minpathsize), 
      "pathways remaining.\n")
  pathway_data = pathway_data[plengths >= minpathsize]
  
  #now run the inner functions depending on method
  if(method=="fc") {
    if(mc.cores > 1){
      #split fcmat into mc.cores matrices and run them in parallel
      splitseq = sort(rep(1:mc.cores, length.out = nrow(FCMAT2))) 
      fclist = split.data.frame(FCMAT2,splitseq, drop = F)
      
      cl = makePSOCKcluster(mc.cores)
      clusterExport(cl, c("pathwayScoreCoreFC"))
      pscorelist = parLapply(cl = cl, X=fclist, fun=pathwayScoreCoreFC,
                               pathset = pathset,
                               dataset = dataset,
                               chem_dict = CHEM_DICT,
                               pathway_data = pathway_data )
      stopCluster(cl)
      #reform output
      pathscoremat = as.data.frame(rbindlist(pscorelist))
    } else {
      pathscoremat = pathwayScoreCoreFC(fcdata = FCMAT2, pathset = pathset,dataset = dataset, 
                                        chem_dict = CHEM_DICT,pathway_data = pathway_data )
    }
    
    #save
    file <- paste0("output/pathway_score_summary/PATHSCOREMAT_",pathset,"_",dataset,"_",method,".RData")
    save(pathscoremat,file=file)
    
  }
  #call gsva scoring
  if(method=="gsva") {
    pathwayScoreCoreGSVA(sk.list,
                         pathset=pathset,
                         dataset=dataset,
                         fcmat=FCMAT2,
                         chem_dict=CHEM_DICT,
                         pathway_data=pathway_data,
                         mc.cores=mc.cores)
  }
  #call mygsea scoring
  if(method=="mygsea") {
    pathwayScoreCoreMYGSEA(sk.list,
                           method = "mygsea",
                           pathset=pathset,
                           dataset=dataset,
                           fcmat=FCMAT2,
                           chem_dict=CHEM_DICT,
                           pathway_data=pathway_data,
                           mc.cores=mc.cores,
                           normalization = T,
                           useranks = T)
  }
  #call mygsea scoring without ranks or normalization
  if(method=="mygsea_norank") {
    pathwayScoreCoreMYGSEA(sk.list,
                           method = "mygsea_norank",
                           pathset=pathset,
                           dataset=dataset,
                           fcmat=FCMAT2,
                           chem_dict=CHEM_DICT,
                           pathway_data=pathway_data,
                           mc.cores=mc.cores,
                           normalization = F,
                           useranks = F)
  }

  #DEPRECATED SCORING FUNCTIONS
  # if(method=="gsvae") {
  #   pathwayScoreCoreGSVAE(sk.list,
  #                        pathset=pathset,
  #                        dataset=dataset,
  #                        fcmat=FCMAT2,
  #                        chem_dict=CHEM_DICT,
  #                        pathway_data=pathway_data,
  #                        mc.cores=mc.cores)
  # }
  # if(method=="gsvap") {
  #   pathwayScoreCoreGSVAPoisson(sk.list,
  #                        pathset=pathset,
  #                        dataset=dataset,
  #                        fcmat=FCMAT2,
  #                        chem_dict=CHEM_DICT,
  #                        pathway_data=pathway_data,
  #                        mc.cores=mc.cores)
  # }
  # if(method=="ssgsea") {
  #   pathwayScoreCoreSSGSEA(sk.list,
  #                        pathset=pathset,
  #                        dataset=dataset,
  #                        fcmat=FCMAT2,
  #                        chem_dict=CHEM_DICT,
  #                        pathway_data=pathway_data,
  #                        mc.cores=mc.cores)
  # }
  # if(method=="plage") {
  #   pathwayScoreCorePLAGE(sk.list,
  #                        pathset=pathset,
  #                        dataset=dataset,
  #                        fcmat=FCMAT2,
  #                        chem_dict=CHEM_DICT,
  #                        pathway_data=pathway_data,
  #                        mc.cores=mc.cores)
  # }
  # if(method=="zscore") {
  #   pathwayScoreCoreZSCORE(sk.list,
  #                         pathset=pathset,
  #                         dataset=dataset,
  #                         fcmat=FCMAT2,
  #                         chem_dict=CHEM_DICT,
  #                         pathway_data=pathway_data,
  #                         mc.cores=mc.cores)
  # }
  # if(method=="fisher") {
  #   pathwayScoreCoreFisher(sk.list,
  #                        pathset=pathset,
  #                        dataset=dataset,
  #                        fcmat=FCMAT2,
  #                        chem_dict=CHEM_DICT,
  #                        pathway_data=pathway_data,
  #                        mc.cores=mc.cores)
  # }
  cat("\nTime taken:",proc.time() - starttime, "\n")

}
