#' Perform the first step of makeing signatures non-redundant
#'
#' @param sigset Name of signature set.
#' @param sigcatalog Nmae of the catalog file
#'
#' @return No output.
#' @export
signatureNonredundant.1 <- function(sigcatalog="signatureDB_master_catalog 2020-05-05",
                                    sigset="screen_large",
                                    cutoff=0.9,
                                    mc.cores=1,
                                    do.load=F) {

  printCurrentFunction()
  if(do.load) {
    cat("load signature data\n")
    file <- paste0("../input/signatures/signatureDB_genelists.RData")
    cat("   ",file,"\n")
    load(file) #genelists
    catalog <- signatureCatalogLoader(sigset,sigcatalog)

    catalog <- catalog[is.element(catalog$signature,names(genelists)),]
    signature_data <- genelists[catalog$signature]
    SDATA <<- signature_data
    GL <<- genelists
  }

  slist <- SDATA
  cat("length: ",length(slist),"\n")
  nlist <- substr(names(slist),1,4)
  mask <- vector(mode="integer",length=length(nlist))
  mask[] <- 1
  mask[nlist=="Rand"] <- 0
  mask[nlist=="CMAP"] <- 0
  slist <- slist[mask==1]
  cat("length: ",length(slist),"\n")

  allsim <- NULL
  for(i in 1:(length(slist))) {
    res <- sigdist.1.all(slist[i],slist[(i+1):length(slist)],mc.cores)
    x <- as.data.frame(unlist(res))
    name.list <- c("s1","s2","similarity")
    y <- as.data.frame(matrix(nrow=nrow(x),ncol=3))
    names(y) <- name.list
    y[,"s1"] <- names(slist)[i]
    y[,"s2"] <- rownames(x)
    y[,"similarity"] <- x[,1]
    y <- y[y$similarity>=cutoff,]
    match <- y[,1]==y[,2]
    y <- y[!match,]
    #cat(names(slist)[i],nrow(y),"\n")
    if(nrow(y)>0) allsim <- rbind(allsim,y)
    if(i%%100==0) cat("finished",i," out of ",length(slist),":",nrow(allsim),"\n")
  }
  file <- paste0("../input/signatures/similarity_",sigset,"_",sigcatalog,".xlsx")
  write.xlsx(allsim,file)
}
sigdist.1.all <- function(s1,slist,mc.cores) {
  #res <- mclapply(X=slist,FUN=sigdist.1.1,s1=s1,mc.cores=mc.cores)
  res <- lapply(X=slist,FUN=sigdist.1.1,s1=s1)
  return(res)
}
sigdist.1.1 <- function(s1,s2) {
  val <- sum(is.element(s1[[1]],s2)) / length(s1[[1]])
  return(val)
}
