#--------------------------------------------------------------------------------------
#' Transpose and filter the fold change matrix FCMAT1 in long format into a
#' gene x sample format.
#'
#' @param chems The CHEMS data frame with chemical information
#' @param method =gene or probe_id
#' @return Global variables are created for the FC matrix (FCMAT2), the SE matrix (SEMAT2)
#' and the chemical dictionary (CHEM_DICT) which translates form the sample key
#' (sample_id_conc_time) to the individual components
#'
#--------------------------------------------------------------------------------------
buildFCMAT2.big <- function(dataset="DMEM_6hr_screen_normal_pe_1",
                            dir="../input/fcdata/",
                            method="gene",do.read=F,
                            dsstox.file="../input/DSSTox/DSSTox_sample_map.xlsx",
                            chemical.file="../input/chemicals/HTTr_screen_sample_map.xlsx") {
  printCurrentFunction(dataset)
  cat("load FCMAT1\n")
  flush.console()
  if(do.read) {
    file <- paste0(dir,"FCMAT1_",dataset,".RData")
    print(file)
    load(file)
    FCMAT1 <<- FCMAT1
    cat("data loaded\n")
  }
  cat("copy FCMAT1 to mat\n")
  flush.console()
  mat <- FCMAT1

  cat("generate the sample key information\n")
  flush.console()
  sid.list <- unique(mat$sample_key)
  temp <- str_split(sid.list,"_")
  sample.map <- as.data.frame(do.call(rbind,temp),stringsAsFactors=F)
  name.list <- c("sample_id","conc_index","conc.string","media","time")
  if(ncol(sample.map)==3) name.list <- c("sample_id","conc_index","conc.string")
  names(sample.map) <- name.list
  conc <- sample.map$conc.string
  units <- conc
  units <- str_replace_all(units,"\\.","")
  for(i in 0:9) units <- str_replace_all(units,as.character(i),"")

  for(unit in unique(units)) conc <- str_replace_all(conc,unit,"")
  conc <- as.numeric(conc)

  sample.map$conc <- conc
  sample.map$units <- units
  sample.map$sample_key <- sid.list
  if(is.element("time",names(sample.map))) {
    time <- sample.map$time
    sample.map$time <- as.numeric(str_replace(time,"h",""))
    sample.map <- sample.map[,c("sample_key","sample_id","conc_index","conc","units","media","time")]
  }
  else sample.map <- sample.map[,c("sample_key","sample_id","conc_index","conc","units")]

  cat("add dtxsid, name, casrn\n")
  smat <- read.xlsx(dsstox.file)
  names(smat)[1] <- "sample_id"
  smat <- unique(smat[,c("sample_id","dtxsid","casrn","name")])
  smat <- smat[!is.na(smat$sample_id),]
  rownames(smat) <- smat$sample_id
  spid.list <- sample.map$sample_id

  cat("check for missing sample_ids\n")
  missing <- spid.list[!is.element(spid.list,smat$sample_id)]
  if(length(missing)>0) {
    cat("missing sample IDs\n")
    print(missing)
    browser()
  }
  cat("write the chemical tables\n")
  temp <- smat[sample.map$sample_id,c("dtxsid","casrn","name")]
  chems <- cbind(sample.map,temp)
  write.xlsx(chems,chemical.file)

  CHEM_DICT <- chems
  file <- paste0(dir,"CHEM_DICT_",dataset,".RData")
  save(CHEM_DICT,file=file)

  cat("build FCMAT2\n")
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
  cat("Dimension of mat3 before filtering: ",dim(mat3),"\n")
  temp <- mat3
  temp[!is.na(mat3)] <- 1
  temp[is.na(mat3)] <- 0
  rs <- rowSums(temp)
  nsample <- dim(mat3)[2]
  mask <- rs
  mask[] <- 1
  mask[rs<0.95*nsample] <- 0
  matfc <- matfc[mask==1,]
  cat("Dimension of mat3 after filtering: ",dim(matfc),"\n")
  matfc <- t(matfc)

  cat("export FCMAT2\n")
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
