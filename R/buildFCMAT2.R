#--------------------------------------------------------------------------------------
#' Transpose and filter the fold change matrix FCMAT1 in long format into a
#' gene x sample format.
#' Function returns FCMAT2 and CHEM_DICT as a list R object
#' @param FCMAT1 The input data structure that contains the data frame from the httrlib buildFCMAT1() function or provided by the user in the same format.
#' FCMAT1 should be a long format table with rows denoting genes/probes and an 'l2fc' column which contains log2 fold-change data from DESeq2 or another method.
#' @param time The time in hours that the chemical dosing was run
#' @param media The name of the media used
#' @param output_dir The directory from which to read all of the raw files
#' @param method Either "gene" or "probe_id". Specifying "gene" will average probes that target the same gene, whereas "probe_id" will keep all probes separate when making the FCMAT2 table.
#' @param sample_type, either "test sample" or "reference chemical" will determine what sample type the FCMAT2/CHEM_DICT tables will include
#' @param writing_to_disk dictates whether all results should be saved to disk (path is defined by output_dir parameter)
#'
#' @return Global variables are created for the FC matrix (FCMAT2)
#' and the chemical dictionary (CHEM_DICT) which translates from the sample key
#' (sample_id_conc_time) to the individual components
#' @importFrom utils flush.console
#' @importFrom reshape2 acast
#' @importFrom stringr str_replace
#' @import data.table
#'
#' @export buildFCMAT2
#' @return list consisting of matrices FCMAT2 and CHEM_DICT
#--------------------------------------------------------------------------------------
buildFCMAT2 <- function(FCMAT1=NULL,
                        time = 24,
                        media = "DMEM",
                        output_dir="../input/fcdata/",
                        method="gene",
                        sample_type="test sample",
                        writing_to_disk=FALSE) {

  printCurrentFunction()

  # if (sample_type=="reference chemical"){
  #   if(do.prep) {
  #     cat("  copy FCMAT1 to mat\n")
  #     mat <- as.data.table(FCMAT1)
  #     sids <- unique(mat$chem_id)
  #     sids.include <- sids[nchar(sids)<=7]
  #     mat <- mat[is.element(mat$chem_id,sids.include),]
  #     MAT2 <<- mat
  #   }
  #   mat <- MAT2
  #   cat("  set up the new sample key: ",nrow(mat),"\n")
  #   sk <- mat$sample_key
  #   sid <- mat$chem_id
  #   pg <- mat$pg_id
  #   sid2 <- paste0(sid,".",pg)
  #   sk2 <- stringr::str_replace(sk,sid,sid2)
  #   mat$chem_id <- sid2
  #   mat$sample_key <- sk2
  # }


  cat("  copy FCMAT1 to mat",nrow(FCMAT1),"\n")
  flush.console()
  mat <- as.data.table(FCMAT1)

  cat("  generate the sample key information\n")
  flush.console()

    if ("casrn" %in% colnames(mat)){
      sample_map <- unique(mat[,c("sample_key","stype", "chem_id","dose_level","conc","conc_unit","dtxsid","chem_name")])
      names(sample_map) <- c("sample_key", "stype", "sample_id","conc_index","conc","units","dtxsid","name")
      sample_map$media <- media
      sample_map$time <- time
  }
  else{
    sample_map <- unique(mat[,c("sample_key", "stype", "chem_id","dose_level","conc","conc_unit","dtxsid","chem_name")])
    names(sample_map) <- c("sample_key", "stype", "sample_id","conc_index","conc","units","dtxsid","name")
    sample_map$media <- media
    sample_map$time <- time
    sample_map$casrn <- "TBD"
  }

  if (sample_type=="test sample"){
    sids <- unique(sample_map[stype == sample_type, sample_id])
    sample_map <- sample_map[sample_id %in% sids,]
    mat <- mat[is.element(mat$sample_key,sample_map$sample_key),]
  }
  else if (sample_type == "reference chemical"){
    sids <- unique(sample_map[stype == sample_type, sample_id])
    sample_map <- sample_map[sample_id %in% sids,]
    mat <- mat[is.element(mat$sample_key,sample_map$sample_key),]
    sk <- mat$sample_key
    sid <- mat$chem_id
    pg <- mat$pg_id
    sid2 <- paste0(sid,".",pg)
    sk2 <- stringr::str_replace(sk,sid,sid2)
    mat$chem_id <- sid2
    mat$sample_key <- sk2

    sample_map <- unique(mat[,c("sample_key","chem_id","dose_level","conc","conc_unit","dtxsid","chem_name")])
    names(sample_map) <- c("sample_key","sample_id","conc_index","conc","units","dtxsid","name")
    sample_map$casrn <- "TBD"
    sample_map$media <- media
    sample_map$time <- time
  }

  CHEM_DICT <- sample_map
  if (writing_to_disk==TRUE){
    file <- paste0(output_dir,"CHEM_DICT_", stringr::str_replace(sample_type, " ", "_"),".RDS")
    saveRDS(CHEM_DICT,file=file)
  }

  if(is.element("gene_symbol",names(mat))) {
    x = names(mat)
    x[is.element(x,"gene_symbol")] = "gene"
    names(mat) = x
  }
  #browser()
  cat("  build FCMAT2\n")
  if(method=="gene") {
    mat2 <- mat[,c("sample_key","gene","l2fc")]
    matfc <- reshape2::acast(mat2,gene~sample_key,mean)

  }
  else if(method=="probe_id") {
    mat2 <- mat[,c("sample_key","probe_id","l2fc")]
    matfc <- reshape2::acast(mat2,probe_id~sample_key,mean)

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
  if (writing_to_disk==TRUE){
    file <- paste0(output_dir,"FCMAT2_", stringr::str_replace(sample_type, " ", "_"),".RDS")
    saveRDS(FCMAT2, file)
  }
  print(dim(FCMAT2))


  flush.console()
  return(list("CHEM_DICT" = CHEM_DICT, "FCMAT2" = FCMAT2))
}
