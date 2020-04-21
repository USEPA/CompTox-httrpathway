#--------------------------------------------------------------------------------------
#'
#' Concatenate the input DESeq2 files
#--------------------------------------------------------------------------------------
concatDESeq2Files <- function(dataset="DMEM_6hr_screen_normal_pe_1",
                              indir="../input/httr_mcf7_screen/meanncnt0_5-plateteffect_0-shrinkage_normal_DMEM_6/",
                              outdir="../input/httr_mcf7_screen/") {
  printCurrentFunction()

  ofile <- paste0(outdir,dataset,".csv")
  file.list <- list.files(indir)
  #file.list <- file.list[1:10]
  nfile <- length(file.list)

  for(i in 1:nfile) {
    file <- file.list[i]
    fname <- paste0(indir,file)
    cat(i," out of ",nfile," : ",file,"\n")
    #print(fname)
    mat <- read.table(fname,header=T,stringsAsFactors=F,sep="\t")
    if(ncol(mat)!=8) {
      cat(ncol(mat),nrow(mat),fname,"\n")
      browser()
    }
    names(mat) <- c("probe_id","sample_key","basemean","l2fc","se","stat","pvalue","padj")
    probes <- as.data.frame(mat[,1])
    names(probes) <- "probe_id"
    genes <- separate(probes,"probe_id",sep="_",into=c("genes","id"))
    mat$gene <- genes[,1]
    mat <- mat[,c("sample_key","probe_id","gene","l2fc","se","basemean","stat","pvalue","padj")]
    if(i==1) write.table(mat,ofile,append=F,row.names=F,sep=",")
    else write.table(mat,ofile,append=T,row.names=F,col.names=F,sep=",")
  }
}
