#--------------------------------------------------------------------------------------
#' Transpose and filter the fold change matrix FCMAT1 in long format into a
#' gene x sample format.
#'
#' @param dataset The name to give to the data set
#' @param dir The directory from which to read all of the raw filesatalog file
#' @param mc.cores The number of cores to use in reading the tsv files
#' @param method Either "gene" or "probe"
#' @param do.read If TRUE, read in the FCMAT1 file and place in a global.
#' @param chemical.file The required map from sample keys to chemical information
#' @param dsstox.file The information mapping chemicals to DSSTox IDs
#/
#' @return Global variables are created for the FC matrix (FCMAT2), the SE matrix (SEMAT2)
#' and the chemical dictionary (CHEM_DICT) which translates form the sample key
#' (sample_id_conc_time) to the individual components
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' u2os_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123_allPG
#' mcf7_ph1_pe1_normal_block_123_excludePG
#' u2os_pilot_pe1_normal_null_full
#' u2os_pilot_pe1_normal_null_pilot
#' u2os_pilot_pe1_normal_null_pilot_lowconc
#' u2os_pilot_pe1_normal_null_pilot_lowconc_lowchem
#'
#--------------------------------------------------------------------------------------
buildFCMAT2.fromDB <- function(dataset="mcf7_ph1_pe1_normal_block_123_excludePG",
                               time=24,
                               media="DMEM",
                               dir="../input/fcdata/",
                               method="gene",
                               do.read=T) {
  printCurrentFunction(dataset)
  cat("  load FCMAT1\n")
  flush.console()
  if(do.read) {
    file <- paste0(dir,"FCMAT1_",dataset,".RData")
    print(file)
    load(file)
    FCMAT1 <<- FCMAT1
    cat("  data loaded\n")
  }
  cat("  copy FCMAT1 to mat",nrow(FCMAT1),"\n")
  flush.console()
  mat <- FCMAT1

  cat("  generate the sample key information\n")
  flush.console()
  sample_map <- unique(mat[,c("sample_key","chem_id","dose_level","conc","conc_unit","dtxsid","chem_name")])
  names(sample_map) <- c("sample_key","sample_id","conc_index","conc","units","dtxsid","name")
  sample_map$casrn <- "TBD"
  sample_map$media <- media
  sample_map$time <- time

  sids <- unique(sample_map$sample_id)
  sids.exclude <- sids[nchar(sids)<=7]
  sample_map <- sample_map[!is.element(sample_map$sample_id,sids.exclude),]
  mat <- mat[is.element(mat$sample_key,sample_map$sample_key),]
  CHEM_DICT <- sample_map
  file <- paste0(dir,"CHEM_DICT_",dataset,".RData")
  save(CHEM_DICT,file=file)

  if(is.element("gene_symbol",names(mat))) {
    x = names(mat)
    x[is.element(x,"gene_symbol")] = "gene"
    names(mat) = x
  }
  #browser()
  cat("  build FCMAT2\n")
  if(method=="gene") {
    mat2 <- mat[,c("sample_key","gene","l2fc")]
    matfc <- acast(mat2,gene~sample_key,mean)

    matp <- mat[,c("sample_key","gene","padj")]
    matpv <- acast(matp,gene~sample_key,mean)
    matpv <- t(matpv)

    mats <- mat[,c("sample_key","gene","se")]
    matse <- acast(mats,gene~sample_key,mean)
    matse <- t(matse)
  }
  else if(method=="probe_id") {
    mat2 <- mat[,c("sample_key","probe_id","l2fc")]
    matfc <- acast(mat2,probe_id~sample_key,mean)

    matp <- mat[,c("sample_key","probe_id","padj")]
    matpv <- acast(matp,probe_id~sample_key,mean)
    matpv <- t(matpv)

    mats <- mat[,c("sample_key","probe_id","se")]
    matse <- acast(mats,probe_id~sample_key,mean)
    matse <- t(matse)

  }

  mat3 <- matfc
  mat3[is.nan(mat3)] <- NA
  cat("  dimension of mat3 before filtering: ",dim(mat3),"\n")
  temp <- mat3
  temp[!is.na(mat3)] <- 1
  temp[is.na(mat3)] <- 0
  rs <- rowSums(temp)
  nsample <- dim(mat3)[2]
  mask <- rs
  mask[] <- 1
  mask[rs<0.95*nsample] <- 0
  matfc <- matfc[mask==1,]
  cat("  dimension of mat3 after filtering: ",dim(matfc),"\n")
  matfc <- t(matfc)

  cat("  export FCMAT2\n")
  FCMAT2 <- matfc
  file <- paste0(dir,"FCMAT2_",dataset,".RData")
  save(FCMAT2,file=file)
  print(dim(FCMAT2))

  matpv <- matpv[,colnames(FCMAT2)]
  matse <- matse[,colnames(FCMAT2)]
  matse.inv <- 1/matse
  matfcse <- matse.inv * FCMAT2

  FCMAT2.FCoverSE <- matfcse
  print(dim(FCMAT2.FCoverSE))
  file <- paste0(dir,"FCMAT2.FCoverSE.",dataset,".RData")
  save(FCMAT2.FCoverSE,file=file)

  FCMAT2.PV <- matpv
  print(dim(FCMAT2.PV))
  file <- paste0(dir,"FCMAT2.PV.",dataset,".RData")
  save(FCMAT2.PV,file=file)

  FCMAT2.SE <- matse
  print(dim(FCMAT2.SE))
  file <- paste0(dir,"FCMAT2.SE.",dataset,".RData")
  save(FCMAT2.SE,file=file)
  flush.console()
}
