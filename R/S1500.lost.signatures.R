#--------------------------------------------------------------------------------------
#' Find the list of signatures that get lost when going to S1500
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#'
#--------------------------------------------------------------------------------------
S1500.lost.signatures <- function(to.file=F,
                                  do.load=F,
                                  dataset="heparg2d_toxcast_pfas_pe1_normal",
                                  sigset="screen_large",
                                  method="fc",
                                  celltype="HepaRG") {
  printCurrentFunction()

  file = "../input/signatures/signatureDB_master_catalog 2021-03-05.xlsx"
  catalog = read.xlsx(file)
  catalog = catalog[catalog[,sigset]==1,]

  file = "../input/S1500/S1500_genes.xlsx"
  s1500.genes = read.xlsx(file)[,1]

  file = "../input/signatures/signatureDB_genelists.RData"
  load(file=file)
  genelists = genelists

  if(do.load) {
    file <- paste0("../input/fcdata/FCMAT2_",dataset,".RData")
    load(file=file)
    FCMAT2 <<- FCMAT2

    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    MAT <<- mat

  }
  all.genes = colnames(FCMAT2)

  sig.list = catalog$signature
  nsig = length(sig.list)

  res = catalog[,c("signature","parent","ngene")]
  res$ngene = 0
  res$ngene.s1500 = 0
  for(i in 1:nrow(res)) {
    sig = res[i,"signature"]
    genes = genelists[sig][[1]]
    x = genes[is.element(genes,all.genes)]
    y = genes[is.element(genes,s1500.genes)]
    res[i,"ngene"] = length(x)
    res[i,"ngene.s1500"] = length(y)
  }
  rownames(res) = res$signature
  res$lost = 0
  plist = unique(res$parent)
  for(parent in plist) {
    slist = res[is.element(res$parent,parent),"signature"]
    counts = res[slist,"ngene.s1500"]
    if(min(counts)<10) res[slist,"lost"] = 1
  }
  cat("lost:",sum(res$lost),"\n")

  res2 = res[,c("parent","lost")]
  res2 = unique(res2[res2$lost==1,])
  slist = res2$parent
  mat = MAT
  mat = mat[mat$top_over_cutoff>4,]
  mat = mat[mat$bmd<1,]
  mat = mat[is.element(mat$signature,slist),]
  file = "../input/S1500/HepaRG_S1500_lost_signature_hits.xlsx"
  write.xlsx(mat,file)
  browser()
}

