#' Read in signature catalog
#'
#' @param sigset Name of the signature set to use; derived from the sigcatalog
#' @param sigcatalog full path to signature catalog xlsx file; default is repo version
#' @param sigdbgenelist full path to signature DB gene list file; default is repo version
#' @importFrom openxlsx read.xlsx
#' @importFrom stringr str_replace
#' @return the trimmed signature table
#' @export signatureCatalogLoader
signatureCatalogLoader <- function(sigset="wgcna",
                                   sigcatalog="../inst/extdata/signatureDB_master_catalog_2022-05-16.xlsx",
                                   sigdbgenelist="../inst/extdata/signatureDB_genelists.RDS") {

  printCurrentFunction()

  file = sigcatalog
  mat <- read.xlsx(file)

  print(sigdbgenelist)
  genelists <- readRDS(sigdbgenelist)


  slist1 <- names(genelists)
  matnorand <- mat[!is.element(mat$source,"Random"),]
  slist2 <- matnorand$signature
  in2not1 <- slist2[!is.element(slist2,slist1)]
  in1not2 <- slist1[!is.element(slist1,slist2)]

  stop <- F
  if(length(in2not1)>0) {
    warning("Signatures are present in the catalog file, but not in the signature database file!")
    print(in2not1)
    stop <- T
  }

  mat <- mat[mat[,sigset]==1,]
  alist <- mat[is.element(mat$type,"nondirectional"),]
  blist <- mat[is.element(mat$type,"directional"),]

  blist <- blist[order(blist$parent),]
  counter <- 0
  for(parent in blist$parent) {
    temp <- blist[is.element(blist$parent,parent),]
    if(nrow(temp)!=2) {
      warning(paste("missing pair:",parent,"\n"))
      counter <- counter+1
    }
  }
  if(counter>0 || stop) {
    warning("Some signature pairs are missing. See warning message above for signatures that need to be fixed.")
  }
  else return(mat)
}

