#--------------------------------------------------------------------------------------
#' Transpose and filter the fold change matrix FCMAT1 in long format into a
#' gene x sample format.
#' This is the method to use when there are conc-response profiles of refchems
#'
#' @param dataset The name to give to the data set
#' @param time The time in hours that the chemical dosing was run
#' @param media THe name of the media used
#' @param dir The directory from which to read all of the raw files
#' @param method Either "gene" or "probe"
#' @param do.read If TRUE, read in the FCMAT1 file and place in a global.
#'
#' @return Global variables are created for the FC matrix (FCMAT2), the SE matrix (SEMAT2)
#' and the chemical dictionary (CHEM_DICT) which translates form the sample key
#' (sample_id_conc_time) to the individual components
#' @export
#--------------------------------------------------------------------------------------
buildFCMAT2.fromDB.refchems <- function(dataset="heparg2d_toxcast_pfas_pe1_normal_v2",
                                        time=24,
                                        media="DMEM",
                                        dir="../input/fcdata/",
                                        method="gene",
                                        do.read=F,
                                        do.prep=T) {
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
  if(do.prep) {
    cat("  copy FCMAT1 to mat\n")
    mat <- FCMAT1
    sids <- unique(mat$chem_id)
    sids.include <- sids[nchar(sids)<=7]
    mat <- mat[is.element(mat$chem_id,sids.include),]
    MAT2 <<- mat
  }
  mat <- MAT2
  cat("  set up the new sample key: ",nrow(mat),"\n")
  sk <- mat$sample_key
  sid <- mat$chem_id
  pg <- mat$pg_id
  sid2 <- paste0(sid,".",pg)
  sk2 <- str_replace(sk,sid,sid2)
  mat$chem_id <- sid2
  mat$sample_key <- sk2
  #browser()

  cat("  generate the sample map\n")
  sample_map <- unique(mat[,c("sample_key","chem_id","dose_level","conc","conc_unit","dtxsid","chem_name")])
  names(sample_map) <- c("sample_key","sample_id","conc_index","conc","units","dtxsid","name")
  sample_map$casrn <- "TBD"
  sample_map$media <- media
  sample_map$time <- time
#  mat <- mat[is.element(mat$sample_key,sample_map$sample_key),]

  CHEM_DICT <- sample_map
  file <- paste0(dir,"CHEM_DICT_",dataset,"_refchems.RData")
  save(CHEM_DICT,file=file)

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
  file <- paste0(dir,"FCMAT2_",dataset,"_refchems.RData")
  save(FCMAT2,file=file)
  print(dim(FCMAT2))

  matpv <- matpv[,colnames(FCMAT2)]
  matse <- matse[,colnames(FCMAT2)]
  matse.inv <- 1/matse
  matfcse <- matse.inv * FCMAT2

  FCMAT2.FCoverSE <- matfcse
  print(dim(FCMAT2.FCoverSE))
  file <- paste0(dir,"FCMAT2.FCoverSE.",dataset,"_refchems.RData")
  save(FCMAT2.FCoverSE,file=file)

  FCMAT2.PV <- matpv
  print(dim(FCMAT2.PV))
  file <- paste0(dir,"FCMAT2.PV.",dataset,"_refchems.RData")
  save(FCMAT2.PV,file=file)

  FCMAT2.SE <- matse
  print(dim(FCMAT2.SE))
  file <- paste0(dir,"FCMAT2.SE.",dataset,"_refchems.RData")
  save(FCMAT2.SE,file=file)
  flush.console()
}
