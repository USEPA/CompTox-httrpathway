#' Gene Concentration Response
#'
#' Wrapper that performs concentration response modeling for gene/probe log2 fold-change values
#'
#' Uses two lowest concentration of each
#' column to estimate noise cutoff (as opposed to signature CR).
#'
#' @param mc.cores Number of parallel cores to use.
#' @param to.file.path when provided, path of RDS file where results are written to
#' @param pval P-value cutoff between 0 and 1.
#' @param aicc If aicc = T, corrected AIC is used instead of first order
#'   (regular) AIC.
#' @param fitmodels Vector of models names to be used. Default is all of them.
#' @param genefile An optional file (.xlsx) that can be used to filter concentration-response modeling for a subset of genes of interest
#' @param FCMAT2 chem/conc by gene or chem/conc by probe. Uses two lowest concentration of each
#'   column to estimate noise cutoff (as opposed to signature CR).
#' @param CHEM_DICT Dataframe with one row per sample key and seven columns:
#'   sample_key, sample_id, conc, time, casrn, name, dtxsid.
#'
#' @import data.table
#' @importFrom parallel makePSOCKcluster clusterExport parLapply stopCluster
#' @importFrom openxlsx read.xlsx
#' @importFrom reshape2 melt
#' @importFrom tcplfit2 concRespCore
#' @importFrom stats sd
#' @return dataframe of concentration response modeling results
#' @export geneConcResp

geneConcResp <- function(mc.cores=20,
                         to.file.path=NULL,
                         pval = .05,
                         aicc = F,
                         fitmodels = c("cnst", "hill", "poly1", "poly2", "pow", "exp2", "exp3",
                                       "exp4", "exp5"),
                         genefile=NULL,
                         FCMAT2,
                         CHEM_DICT) {

  printCurrentFunction()
  starttime = proc.time()

  if(!is.null(genefile)) {
    temp = read.xlsx(genefile)
    temp = temp[temp[,2]>2,]
    genelist = temp[,1]
    genelist = genelist[is.element(genelist,colnames(FCMAT2))]
    FCMAT2 = FCMAT2[,genelist]
  }


  sk.list = CHEM_DICT$sample_key
  FCMAT2 = FCMAT2[sk.list,]
  rownames(CHEM_DICT) = CHEM_DICT$sample_key
  CHEM_DICT0 = CHEM_DICT
  FCMAT20 = FCMAT2
  ################################################################
  # dlist = unique(CHEM_DICT$dtxsid)
  # CHEM_DICT = CHEM_DICT[is.element(CHEM_DICT$dtxsid,dlist),]
  # nr = nrow(CHEM_DICT)
  # FCMAT2 = FCMAT2[rownames(CHEM_DICT),1:500]
  ################################################################

  cat("files read\n")
  print(dim(FCMAT2))
  #melt FCMAT to create one row per value data table
  cat("reshape\n")
  genemat = melt(FCMAT2, value.name = "l2fc",variable.factor = F)
  genemat[,1] = as.character(genemat[,1])
  genemat[,2] = as.character(genemat[,2])
  colnames(genemat)[1:2] = c("sample_key", "gene")
  genemat$gene = as.character(genemat$gene)
  genemat$sample_key = as.character(genemat$sample_key)
  cat("merge\n")

  genetab = data.table(genemat,key="sample_key")
  chemtab = data.table(CHEM_DICT,key="sample_key")
  temptab = merge(genetab,chemtab)
  temp = as.data.table(temptab)
  genemat = temp
  cat("remove NAs\n")
  genemat = genemat[!is.na(genemat$l2fc),]
  cat("genemat built\n")

  #get noise estimate from two lowest concs (could alternatively use semat)
  lowresps = genemat[genemat$conc <= .1, c("gene", "l2fc")]
  stats = setDT(lowresps)[, list(cutoff = quantile(abs(l2fc),1-pval), onesd = sd(l2fc), bmed = 0), by = list(gene)]
  genemat$cutoff = stats$cutoff[match(genemat$gene, stats$gene)]
  genemat$bmed = 0#stats$bmed[match(genemat$gene, stats$gene)]
  genemat$onesd = stats$onesd[match(genemat$gene, stats$gene)]

  #aggregate genemat by unique sample/signature per row; data table is considerably faster than aggregate
  genemat = setDT(genemat)[, list(conc = list(conc),resp = list(l2fc)),
                           by = list(sample_id, dtxsid, casrn, name, time, gene, bmed, cutoff, onesd)]

  ordering = order(tolower(genemat$name), tolower(genemat$gene))
  genemat = genemat[ordering,]

  genemat = as.list(as.data.frame(t(genemat), stringsAsFactors = F))
  cat("genemat as list\n")

  #loop through signatureConcRespcore_pval, which is generic for any conc/resp
  if(mc.cores > 1){
    cat("start running multi core\n")
    cl = makePSOCKcluster(mc.cores)
    clusterExport(cl, c("acy", "acgnlsobj", "bmdbounds", "bmdobj",  "cnst", "exp2", "exp3", "exp4", "exp5", "fitcnst", "fithill", "fitgnls",
                        "fitcnst", "fitpoly1", "fitpoly2", "fitpow", "fitexp2", "fitexp3","fitexp4", "fitexp5",
                        "gnls" , "gnlsderivobj", "hillfn", "hitcontinner","hitloginner","loggnls", "loghill", "nestselect",
                        "poly1", "poly2", "pow", "tcplObj", "toplikelihood"))
    GENE_CR = parLapplyLB(cl = cl, X=genemat, fun=concRespCore, fitmodels = fitmodels,
                          aicc = aicc, chunk.size = ceiling(length(genemat)/5/mc.cores) )
  } else {
    cat("start running single core\n")
    GENE_CR = lapply(X=genemat, FUN = concRespCore, fitmodels = fitmodels, aicc = aicc)
  }
  cat("finish running\n")

  #reformat and save
  GENE_CR = as.data.frame(rbindlist(GENE_CR))
  rm(genemat)

  if(!is.null(to.file.path)){
    saveRDS(GENE_CR,to.file.path)
  }
  cat("finish writing\n")

  if(mc.cores > 1) stopCluster(cl)
  print(proc.time() - starttime)

  return(GENE_CR)

}
