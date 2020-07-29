#' Perform the second step of makeing signatures non-redundant
#'
#' @param sigset Name of signature set.
#' @param sigcatalog Nmae of the catalog file
#'
#' @return No output.
#' @export
signatureNonredundant.2 <- function(sigcatalog="signatureDB_master_catalog 2020-05-05",
                                    sigset="screen_large",
                                    do.load=F) {

  printCurrentFunction()
  if(do.load) {
    cat("load signature data\n")
    file <- paste0("../input/signatures/signatureDB_genelists.RData")
    cat("   ",file,"\n")
    load(file) #genelists
    catalog <- signatureCatalogLoader(sigset,sigcatalog)

    catalog <- catalog[is.element(catalog$signature,names(genelists)),]
    CATALOG <<- catalog
    signature_data <- genelists[catalog$signature]
    SDATA <<- signature_data
    GL <<- genelists
  }
  catalog <- CATALOG
  rownames(catalog) <- catalog$signature
  file <- paste0("../input/signatures/similarity_",sigset,"_",sigcatalog,".xlsx")
  mat <- read.xlsx(file)

  slist <- c(mat[,1],mat[,2])
  temp <- substr(slist,1,4)
  mask <- vector(mode="integer",length=length(temp))
  mask[] <- 1
  mask[temp=="Rand"] <- 0
  slist <- slist[mask==1]
  res <- as.data.frame(table(slist))
  names(res) <- c("signature","count")
  res$ngene <- catalog[res$signature,"ngene"]
  res <- res[order(res$count,decreasing=T),]
  file <- paste0("../input/signatures/similarity_",sigset,"_",sigcatalog," step 2.xlsx")
  write.xlsx(res,file)
}
