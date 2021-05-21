#' Gene Concentration Response
#'
#' Wrapper that performs concentration response modeling for gene or probe l2fc's
#'
#' If conthits = T and nametag is NULL, nametag will be set to "conthits". Loads
#' an FCMAT2 and CHEM_DICT corresponding to given dataset. FCMAT should be
#' chem/conc by gene or chem/conc by probe. Uses two lowest concentration of each
#' column to estimate noise cutoff (as opposed to signature CR). Also, doesn't
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
#'
#'  heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123_allPG
#'  u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#' u2os_pilot_pe1_normal_null_pilot_lowconc
#'
#' DMEM_6hr_pilot_normal_pe_1 - MCF7 pilot
#'
#'  u2os_toxcast_pfas_pe1_normal_refchems
#' heparg2d_toxcast_pfas_pe1_normal_refchems
#'

geneConcResp <- function(dataset="heparg2d_toxcast_pfas_pe1_normal_refchems",
                         mc.cores=20,
                         to.file=T,
                         pval = .05,
                         nametag = NULL,
                         conthits = T,
                         aicc = F,
                         fitmodels = c("cnst", "hill", "poly1", "poly2", "pow", "exp2", "exp3",
                                       "exp4", "exp5"),
                         genefile=NULL
                         ) {

  printCurrentFunction(dataset)
  starttime = proc.time()
  if(is.null(nametag) && conthits) nametag = "conthits"
  if(!is.null(nametag)) nametag = paste0("_", nametag)

  #get FCMAT and CHEM_DICT
  file <- paste0("../input/fcdata/FCMAT2_",dataset,".RData")
  load(file)
  file <- paste0("../input/fcdata/CHEM_DICT_",dataset,".RData")
  load(file)

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
  genemat = reshape2::melt(FCMAT2, value.name = "l2fc",variable.factor = F)
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
                          conthits =conthits, aicc = aicc, chunk.size = ceiling(length(genemat)/5/mc.cores) )
  } else {
    cat("start running single core\n")

    # GENE_CR = NULL
    # for(i in 1:length(genemat)) {
    #   resp = concRespCore(genemat[[i]],fitmodels,force.fit=T,verbose=T,aicc=aicc)
    #   GENE_CR = rbind(GENE_CR,resp)
    #   browser()
    # }

    GENE_CR = lapply(X=genemat, FUN = concRespCore, fitmodels = fitmodels, conthits= conthits, aicc = aicc)
  }
  cat("finish running\n")

  #reformat and save
  GENE_CR = as.data.frame(rbindlist(GENE_CR))
  rm(genemat)

  if(to.file){
    file <- paste0("../output/gene_conc_resp_summary/GENE_CR_",dataset,"_", pval, nametag ,".RData")
    save(GENE_CR,file=file)
  } else return(GENE_CR)
  cat("finish writing\n")

  if(mc.cores > 1) stopCluster(cl)
  print(proc.time() - starttime)

}
