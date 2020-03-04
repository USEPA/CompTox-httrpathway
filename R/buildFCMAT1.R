#--------------------------------------------------------------------------------------
#' Build the FCMAT1 data set
#'
#' @param dataset The name to give to the data set
#' @param dir The directory from which to read all of the raw filesatalog file
#' @param filetype Either tsv or RData
#' @return A file with the FCMAT1 data is written to "../input/fcdata/FCMAT0_",dataset,".RData"
#'
#'
#' buildFCMAT1(dataset="DMEM_6hr_screen_normal_pe_1",dir="../input/httr_mcf7_screen/meanncnt0_5-plateteffect_0-shrinkage_normal_DMEM_6/")
#'
#' http://www.win-vector.com/blog/2015/07/efficient-accumulation-in-r/
#'
#' @export
#--------------------------------------------------------------------------------------
buildFCMAT1 = function(dataset="DMEM_6hr_pilot_normal_pe_1",
                       dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_0-shrinkage_normal/DMEM_6/",
                       filetype="tsv",
                       version=1){
  printCurrentFunction()

  if(version==1) {
    file.list <- list.files(dir,full.names=TRUE)
    #file.list <- file.list[1:10]
    All <- lapply(file.list,function(i){
      read.table(i,header=T,stringsAsFactors=F,sep="\t")
    })
    mat <- do.call(rbind.data.frame, All)
    names(mat) <- c("probe_id","sample_key","basemean","l2fc","se","stat","pvalue","padj")
    probes <- as.data.frame(mat[,1])
    names(probes) <- "probe_id"
    genes <- separate(probes,"probe_id",sep="_",into=c("genes","id"))
    mat$gene <- genes[,1]
    omat <- mat
    #browser()
  }

  if(version==0) {
    file.list <- list.files(dir)
    omat <- NULL
    counter <- 0
    count <- length(file.list)
    for(file in file.list) {
      fname <- paste0(dir,file)
      counter <- counter+1
      cat(counter," out of ",count," : ",file,"\n")
      print(fname)
      if(filetype=="tsv") {
        mat <- read.table(fname,header=T,stringsAsFactors=F,sep="\t")
      }
      else if(filetype=="RData") {
        load(fname)
        mat <- FCests
      }
      else {
        cat("invalid file type",filetype,"\n")
        browser()
      }

      names(mat) <- c("probe_id","sample_key","basemean","l2fc","se","stat","pvalue","padj")
      probes <- as.data.frame(mat[,1])
      names(probes) <- "probe_id"
      genes <- separate(probes,"probe_id",sep="_",into="genes")
      mat$gene <- genes[,1]
      omat <- rbind(omat,mat)
    }
  }
  FCMAT1 <- omat
  file <- paste0("../input/fcdata/FCMAT1_",dataset,".RData")
  save(FCMAT1,file=file)
}
