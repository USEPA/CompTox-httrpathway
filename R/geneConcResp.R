#' Gene Concentration Response
#'
#' Wrapper that performs concentration response modeling for gene/probe log2 fold-change values
#'
#' Uses two lowest concentration of each
#' column to estimate noise cutoff (as opposed to signature CR).
#'
#' @param mc.cores Number of parallel cores to use.
#' @param to.file.path when provided, path of RDS file where results are written to
#' @param pval Desired cutoff p-value.
#' @param nlowconc Only include the lowest nlowconc concentrations for each chemical when estimating the noise band
#' @param aicc aicc = T uses corrected AIC to choose winning method; otherwise regular AIC
#' @param fitmodels Vector of models names to be used. Default is "cnst", "hill", "poly1", "poly2",
#' "pow", "exp2", "exp3", "exp4", "exp5"
#' @param genefile An optional file (.xlsx) that can be used to filter concentration-response modeling for a subset of genes of interest
#' @param FCMAT2 Sample by gene matrix of log2(fold change)'s. Rownames are
#'   sample keys and colnames are genes.
#' @param CHEM_DICT Dataframe with one row per sample key and seven columns:
#'   sample_key, sample_id, conc, time, casrn, name, dtxsid.
#'
#' @import data.table
#' @import tcplfit2
#' @importFrom parallel makePSOCKcluster clusterExport parLapply stopCluster
#' @importFrom openxlsx read.xlsx
#' @importFrom reshape2 melt
#' @importFrom stats sd
#' @return dataframe of concentration response modeling results
#' @export geneConcResp

geneConcResp <- function(mc.cores=20,
                         to.file.path=NULL,
                         pval = .05,
                         nlowconc = 2,
                         aicc = F,
                         fitmodels = c("cnst", "hill", "poly1", "poly2", "pow", "exp2", "exp3",
                                       "exp4", "exp5"),
                         genefile=NULL,
                         FCMAT2,
                         CHEM_DICT) {

  printCurrentFunction()
  starttime = proc.time()

  if(!is.null(genefile)) {
    temp <- read.xlsx(genefile)
    temp <- temp[temp[,2]>2,]
    genelist <- temp[,1]
    genelist <- genelist[is.element(genelist,colnames(FCMAT2))]
    FCMAT2 <- FCMAT2[,genelist]
  }


  sk.list <- CHEM_DICT$sample_key
  FCMAT2 <- FCMAT2[sk.list,]
  rownames(CHEM_DICT) <- CHEM_DICT$sample_key

  ################################################################
  # dlist <- unique(CHEM_DICT$dtxsid)
  # CHEM_DICT <- CHEM_DICT[is.element(CHEM_DICT$dtxsid,dlist),]
  # nr <- nrow(CHEM_DICT)
  # FCMAT2 <- FCMAT2[rownames(CHEM_DICT),1:500]
  ################################################################

  cat("files read\n")
  print(dim(FCMAT2))
  #melt FCMAT to create one row per value data table
  cat("reshape\n")
  genemat <- reshape2::melt(FCMAT2, value.name = "l2fc",variable.factor = F)
  genemat[,1] <- as.character(genemat[,1])
  genemat[,2] <- as.character(genemat[,2])
  colnames(genemat)[1:2] = c("sample_key", "gene")
  genemat$gene <- as.character(genemat$gene)
  genemat$sample_key <- as.character(genemat$sample_key)
  cat("merge\n")

  genetab <- data.table(genemat,key="sample_key")
  chemtab <- data.table(CHEM_DICT,key="sample_key")
  genetab <- merge(genetab,chemtab)
  genetab <- as.data.table(genetab)
  genetab
  cat("remove NAs\n")
  genetab <- genetab[!is.na(genetab$l2fc),]
  cat("genetab built\n")

  #get noise estimate from two lowest concs (could alternatively use semat)
  lowresps <- genetab[conc_index <= nlowconc, c("gene", "l2fc")]
  stats <- setDT(lowresps)[, list(cutoff = quantile(abs(l2fc),1-pval), onesd = sd(l2fc), bmed = 0), by = list(gene)]
  genetab$cutoff <- stats$cutoff[match(genetab$gene, stats$gene)]
  genetab$bmed <- 0#stats$bmed[match(genetab$gene, stats$gene)]
  genetab$onesd <- stats$onesd[match(genetab$gene, stats$gene)]

  #aggregate genetab by unique sample/signature per row; data table is considerably faster than aggregate
  genetab <- setDT(genetab)[, list(conc = list(conc),resp = list(l2fc)),
                            by = list(sample_id, dtxsid, casrn, name, time, gene, bmed, cutoff, onesd)]

  ordering <- order(tolower(genetab$name), tolower(genetab$gene))
  genetab <- genetab[ordering,]

  genetab <- as.list(as.data.frame(t(genetab), stringsAsFactors = F))
  cat("genetab as list\n")

  #loop through signatureConcRespcore_pval, which is generic for any conc/resp
  if(mc.cores > 1){
    cat("start running multi core\n")
    cl <- makePSOCKcluster(mc.cores)
    clusterExport(cl, c("acy", "acgnlsobj", "bmdbounds", "bmdobj",  "cnst", "exp2", "exp3", "exp4", "exp5", "fitcnst", "fithill", "fitgnls",
                        "fitcnst", "fitpoly1", "fitpoly2", "fitpow", "fitexp2", "fitexp3","fitexp4", "fitexp5",
                        "gnls" , "gnlsderivobj", "hillfn", "hitcontinner","hitloginner","loggnls", "loghill", "nestselect",
                        "poly1", "poly2", "pow", "tcplObj", "toplikelihood"))
    GENE_CR <- parLapplyLB(cl = cl, X=genetab, fun=tcplfit2::concRespCore, fitmodels = fitmodels,
                           aicc = aicc, chunk.size = ceiling(length(genetab)/5/mc.cores) )
  } else {
    cat("start running single core\n")
    GENE_CR <- lapply(X=genetab, FUN = tcplfit2::concRespCore, fitmodels = fitmodels, aicc = aicc)
  }
  cat("finish running\n")

  #reformat and save
  GENE_CR <- as.data.frame(rbindlist(GENE_CR))
  rm(genetab)

  if(!is.null(to.file.path)){
    saveRDS(GENE_CR,to.file.path)
  }
  cat("finish writing\n")

  if(mc.cores > 1) stopCluster(cl)
  print(proc.time() - starttime)

  return(GENE_CR)

}
