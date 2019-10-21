#' Gene Concentration Response
#' 
#' Wrapper that performs concentration response modeling for gene or probe l2fc's
#' 
#' If conthits = T and nametag is NULL, nametag will be set to "conthits". Loads
#' an FCMAT2 and CHEM_DICT corresponding to given dataset. FCMAT should be 
#' chem/conc by gene or chem/conc by probe. Uses two lowest concentration of each
#' column to estimate noise cutoff (as opposed to pathway CR). Also, doesn't
#' currently contain a plotting option.
#'
#' @param dataset String that identifies data set.
#' @param mc.cores Number of parallel cores to use.
#' @param to.file If TRUE, results are written to an RData file, otherwise they
#'   are returned.
#' @param pval P-value cutoff between 0 and 1.
#' @param nametag Optional identifier attached to the output name that usually
#'   is used to signify that an unusual option was used.
#' @param conthits If conthits = T, continuous hitcalls are calculated; otherwise
#'   discrete hitcalls are used.
#' @param aicc If aicc = T, corrected AIC is used insstead of first order 
#'   (regular) AIC.
#' @param fitmodels Vector of models names to be used. Default is all of them. 
#' 
#' @import data.table
#' @import parallel
#'
#' @return If to.file = F, data frame containing results; otherwise, nothing.
#' @export
geneConcResp <- function(dataset="ph1_100normal_pid",
                                 mc.cores=1,
                                 to.file=F,
                                 pval = .05,
                                 nametag = NULL,
                                 conthits = F,
                                 aicc = F,
                                 fitmodels = c("cnst", "hill", "gnls", "poly1", "poly2", "pow", "exp2", "exp3", 
                                               "exp4", "exp5")) {
  
  starttime = proc.time()

  if(is.null(nametag) && conthits) nametag = "conthits"
  if(!is.null(nametag)) nametag = paste0("_", nametag)
  
  #get FCMAT and CHEM_DICT
  file <- paste0("input/fcdata/FCMAT2_",dataset,".RData")
  load(file)
  
  file <- paste0("input/fcdata/CHEM_DICT_",dataset,".RData")
  load(file)
  
  #melt FCMAT to create one row per value data table
  genemat = melt(FCMAT2, value.name = "l2fc",variable.factor = F)
  colnames(genemat)[1:2] = c("sample_key", "gene")
  genemat = genemat[!is.na(genemat$l2fc),]
  genemat = cbind(genemat, CHEM_DICT[genemat$sample_key,colnames(CHEM_DICT) != "sample_key"])
  
  #get noise estimate from two lowest concs (could alternatively use semat)
  lowresps = genemat[genemat$conc <= .1, c("gene", "l2fc")]
  stats = setDT(lowresps)[, list(cutoff = quantile(abs(l2fc),1-pval), onesd = sd(l2fc), bmed = median(l2fc)), by = list(gene)]
  genemat$cutoff = stats$cutoff[match(genemat$gene, stats$gene)]
  genemat$bmed = stats$bmed[match(genemat$gene, stats$gene)]
  genemat$onesd = stats$onesd[match(genemat$gene, stats$gene)]
  
  #aggregate genemat by unique sample/pathway per row; data table is considerably faster than aggregate
  genemat = setDT(genemat)[, list(conc = list(conc),resp = list(l2fc)), 
                                     by = list(sample_id, dsstox_substance_id, casrn, name, time, gene, bmed, cutoff, onesd)]
  
  ordering = order(tolower(genemat$name), tolower(genemat$gene))
  genemat = genemat[ordering,]
  
  genemat = as.list(as.data.frame(t(genemat), stringsAsFactors = F))
  
  #loop through pathwayConcRespcore_pval, which is generic for any conc/resp
  if(mc.cores > 1){
    cl = makePSOCKcluster(mc.cores)
    clusterExport(cl, c("acy", "acgnlsobj", "bmdbounds", "bmdobj", "cnst", "exp2", "exp3", "exp4", "exp5", "fitcnst", "fithill", "fitgnls", 
                        "fitcnst", "fitpoly1", "fitpoly2", "fitpow", "fitexp2", "fitexp3","fitexp4", "fitexp5", 
                        "gnls" , "gnlsderivobj", "hillfn", "hitcontinner","hitloginner","httrFit", "loggnls", "loghill", "nestselect",
                        "poly1", "poly2", "pow", "tcplObj", "toplikelihood"))
    GENE_CR = parLapplyLB(cl = cl, X=genemat, fun=pathwayConcRespCore_pval, fitmodels = fitmodels, 
                             conthits =conthits, aicc = aicc, chunk.size = ceiling(length(genemat)/5/mc.cores) )
  } else {
    GENE_CR = lapply(X=genemat, FUN = pathwayConcRespCore_pval, fitmodels = fitmodels, conthits= conthits, aicc = aicc)
  }
  
  #reformat and save
  GENE_CR = as.data.frame(rbindlist(GENE_CR))
  rm(genemat)
  
  if(to.file){
    file <- paste0("output/gene_conc_resp_summary/GENE_CR_",dataset,"_", pval, nametag ,".RData")
    save(GENE_CR,file=file)
  } else return(GENE_CR)
  
  if(mc.cores > 1) stopCluster(cl)
  print(proc.time() - starttime)

}