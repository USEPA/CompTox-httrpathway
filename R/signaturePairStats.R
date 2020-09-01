#--------------------------------------------------------------------------------------
#' find the overlap between pairs of signaturss that occur often
#'
#' @param dataset The L2fc matrix data set
#' @param dir The directory where the data file lives
#' @param do.read If TRUE, read in FCMAT2 to a gloabal
#' @return
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#' u2os_toxcast_pfas_pe1_normal
#--------------------------------------------------------------------------------------
signaturePairStats <- function(do.read=F,
                               dataset="heparg2d_toxcast_pfas_pe1_normal",
                               sigset="screen_large",
                               method="fc",
                               celltype="HepaRG",
                               sigcatalog="signatureDB_master_catalog 2020-08-14") {
  printCurrentFunction(paste(dataset,sigset,method))
  if(do.read) {
    file <- paste0("../input/signatures/",sigcatalog,".xlsx")
    catalog <- read.xlsx(file)
    catalog <- catalog[catalog[,sigset]==1,]
    CATALOG <<- catalog

    file = paste0("../input/signatures/signatureDB_genelists.RData")
    print(file)
    load(file=file)
    GENELISTS <<- genelists[catalog$signature]

    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat <- SIGNATURE_CR
    mat <- mat[mat$bmd<1,]
    mat <- mat[mat$hitcall>0.99,]
    mat <- mat[mat$top_over_cutoff>3,]
    MAT <<- mat
  }

  mat = MAT
  #mat <- mat[1:10000,]
  sig1 <- NULL
  sig2 <- NULL
  clist <- sort(unique(mat$name))
  for(chem in clist) {
    matc <- mat[is.element(mat$name,chem),]
    matc$auc <- matc$top * (3-log10(matc$bmd))
    matc <- matc[order(matc$auc,decreasing=T),]
    temp <- matc[,"signature"]
    temp <- sort(temp)
    cat(chem,":",length(temp),"\n")
    if(length(temp)>100) temp <- temp[1:100]
    if(length(temp)>1) {
      for(i in 1:(length(temp)-1)) {
        for(j in (i+1):length(temp)) {
          sig1 <- c(sig1,temp[i])
          sig2 <- c(sig2,temp[j])
        }
      }
    }
  }
  sig <- paste0(sig1,"|",sig2)
  tsig <- table(sig)
  ntsig <- length(tsig)
  name.list <- c("sig1","sig2","n1","n2","n12","copies","similarity")
  res <- as.data.frame(matrix(nrow=ntsig,ncol=length(name.list)))
  names(res) <- name.list
  for(i in 1:ntsig) {
    pair <- names(tsig)[i]
    pairs <- str_split(pair,"\\|")[[1]]
    sig1 <- pairs[1]
    sig2 <- pairs[2]
    res[i,"sig1"] <- sig1
    res[i,"sig2"] <- sig2
    g1 <- unique(sort(as.character(unlist(GENELISTS[sig1]))))
    g2 <- unique(sort(as.character(unlist(GENELISTS[sig2]))))
    res[i,"copies"] <- as.integer(tsig[i])
    n1 <- length(g1)
    n2 <- length(g2)
    n12 <- length(g1[is.element(g1,g2)])
    if(n12>n1 || n12>n2) browser()
    sim12 <- max(n12/n1,n12/n2)
    if(is.nan(sim12)) sim12 <- 0
    res[i,"n1"] <- n1
    res[i,"n2"] <-n2
    res[i,"n12"] <- n12
    res[i,"similarity"] <- sim12
  }

  file <- paste0("../output/signature_cluster/",celltype,"/signature_pair_stats_",celltype,"_",sigset,"_",dataset,"_",method,".xlsx")
  print(file)
  write.xlsx(res,file)
}
