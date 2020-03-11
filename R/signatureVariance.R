#--------------------------------------------------------------------------------------
#' Visualize the variance by signature
#'
#' @param min.ngene Signatures will only be saved if the number of genes is >= this value
#' @param max.ngene Signatures will only be saved if the number of genes is <= this value
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
signatureVariance <- function(to.file=F,
                              do.load=F,
                              sigset="pilot_large_all_CMAP",
                              sigcatalog="signatureDB_master_catalog 2020-01-31",
                              basedir="../input/fcdata/",
                              dataset="DMEM_6hr_screen_normal_pe_1",
                              name="Fulvestrant",
                              parent.list=c("CMAP fulvestrant 1e-08 100 985 100","CMAP estradiol 1e-08 100 8350 100"),
                              conc.list=c()){
  printCurrentFunction()

  if(do.load) {
    file <- paste0(basedir,"FCMAT2_",dataset,".RData")
    print(file)
    load(file)
    file <- paste0(basedir,"CHEM_DICT_",dataset,".RData",sep="")
    print(file)
    load(file)
    rownames(CHEM_DICT) <- CHEM_DICT[,"sample_key"]

    FCMAT2 <<- FCMAT2
    CHEM_DICT <<- CHEM_DICT
    cat("   load signature data\n")
    file <- paste0("../input/signatures/signatureDB_genelists.RData")
    cat("   ",file,"\n")
    load(file) #genelists
    genelists <<- genelists

    annotations <<- signatureCatalogLoader(sigset,sigcatalog)
  }

  genelists.use <- genelists[annotations$signature]
  FCMAT2[is.nan(FCMAT2)] <- 0
  ns <- length(genelists.use)
  #ns <- 200
  name.list <- c("signature","ngene","mean","sd","min","max","med","mad")
  res <- as.data.frame(matrix(nrow=length(ns),ncol=length(name.list)))
  names(res) <- name.list

  for(i in 1:ns) {
    sig <- names(genelists.use)[i]
    gl <- genelists.use[i][[1]]
    gl <- gl[is.element(gl,colnames(FCMAT2))]
    temp <- FCMAT2[,gl]
    x <- as.numeric(temp)
    x <- sort(x)
    res[i,"signature"] <- sig
    res[i,"ngene"] <- length(gl)
    res[i,"mean"] <- mean(x)
    res[i,"sd"] <- sd(x)
    res[i,"med"] <- median(x)
    res[i,"mad"] <- mad(x)
    res[i,"min"] <- min(x)
    res[i,"max"] <- max(x)
    if(i%%100==0) cat("finished ",i,"out of ",ns,"\n")
  }

  file <- paste0("../output/signature_corr/signatureVariance.xlsx")
  write.xlsx(res,file)
  if(to.file) {
    fname <- paste0("../output/signature_corr/signatureVariance.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(4,2),mar=c(4,4,2,2))
 # hist(res$ngene,main="ngene")
  plot(res$sd~res$mean,xlab="mean",ylab="sd")
  plot(res$mad~res$med,xlab="med",ylab="mad")
  hist(res$min,main="min")
  hist(res$max,main="max")
  hist(res$mean,main="mean")
  hist(res$sd,main="sd")
  hist(res$med,main="med")
  hist(res$mad,main="mad")

   if(!to.file) browser()
  else dev.off()
}
