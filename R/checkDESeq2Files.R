#--------------------------------------------------------------------------------------
#'
#' Check the DESeq2 input files
#--------------------------------------------------------------------------------------
checkDESeq2Files <- function(dir="../input/httr_mcf7_screen/meanncnt0_5-plateteffect_0-shrinkage_normal_DMEM_6/") {
  printCurrentFunction()
  file.list <- list.files(dir)
  omat <- NULL
  counter <- 0
  count <- length(file.list)
  temp <- NULL
  for(file in file.list) {
    fname <- paste0(dir,file)
    counter <- counter+1
    cat(counter," out of ",count," : ",file,"\n")
    #print(fname)
    mat <- read.table(fname,header=T,stringsAsFactors=F,sep="\t")
    if(ncol(mat)!=8) {
      cat(ncol(mat),nrow(mat),fname,"\n")
      browser()
    }
  }
 }
