#--------------------------------------------------------------------------------------
#' Build the FCMAT1 data set
#' This is a variant on the original buildFCMAT1 to handle large data sets.
#'
#' @param dataset The name to give to the data set
#' @param dir The directory from which to read all of the raw filesatalog file
#' @param filetype Either tsv or RData
#' @return A file with the FCMAT1 data is written to "../input/fcdata/FCMAT1_",dataset,".RData"
#'
#'
#' buildFCMAT1.big(dataset="DMEM_6hr_screen_normal_pe_1",dir="../input/httr_mcf7_screen/meanncnt0_5-plateteffect_0-shrinkage_normal_DMEM_6/")
#'
#' http://www.win-vector.com/blog/2015/07/efficient-accumulation-in-r/
#'
#' @export
#--------------------------------------------------------------------------------------
buildFCMAT1.big = function(dataset="DMEM_6hr_screen_normal_pe_1",
                           dir="../input/httr_mcf7_screen/meanncnt0_5-plateteffect_0-shrinkage_normal_DMEM_6/",
                           use.gz=F){
  printCurrentFunction()


  # concatenate all of the tsv files
  # this creates a large file called ...
  ### concatDESeq2Files(dataset=dataset,indir=dir,outdir="../input/httr_mcf7_screen/")

  ofile <- paste0("../input/httr_mcf7_screen/",dataset,".csv")
  # gzip this file so that it's easier to read in
  # do this on the command line
  # gzip ofile

  # read in the gzipped file and write out as an RData file
  if(use.gz) {
    ifile <- paste0(ofile,".gz")
    cat("start reading zip file\n")
    zz=gzfile(ifile,'rt')
    cat("finish reading zip file and start unzipping\n")
    omat=read.csv(zz,header=F,stringsAsFactors=F)
    omat=read.csv(zz,header=F,stringsAsFactors=F)
    cat("Finish conversion to unzipped version\n")
  }
  else {
    cat("start reading csv file\n")
    omat=read.csv(ofile,header=F,stringsAsFactors=F)
    cat("finish reading csv file\n")
  }
  print(names(omat))
  print(dim(omat))
  name.list <- c("sample_key","probe_id","gene","l2fc","se","basemean","stat","pvalue","padj")
  names(omat) <- name.list
  cat("change the column types\n")
  omat <- omat[2:nrow(omat),]
  for(i in 1:3) omat[,i] <- as.character(omat[,i])
  for(i in 4:9) omat[,i] <- as.numeric(omat[,i])

  FCMAT1 <- omat
  file <- paste0("../input/fcdata/FCMAT1_",dataset,".RData")
  cat("start saving RData file\n")
  save(FCMAT1,file=file)
  cat("finish saving RData file\n")
}
