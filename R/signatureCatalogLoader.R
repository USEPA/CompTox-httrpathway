#' Merge the up and down halves of the pathway data
#'
#' @param sigset Name of the signature set.
#' @param sigcatlog Nmae of the catalog file
#' @return the trimmed signature table
#' @export
signatureCatalogLoader <- function(sigset="wgcna",
                                   sigcatalog="signatureDB_wgcna_mcf7_ph1_pe1_normal_good_pg_MCF7_12_10_catalog") {

  printCurrentFunction()

  file = paste0("../input/signatures/",sigcatalog,".xlsx")
  mat <- read.xlsx(file)
  if(sigset!="wgcna") {
    file = paste0("../input/signatures/signatureDB_genelists.RData")
    print(file)
    load(file=file)
  }
  else {
    name <- str_replace(sigcatalog,"_catalog","")
    file = paste0("../input/signatures/",name,".RData")
    print(file)
    load(file=file)
    genelists <- sigdb
  }
  #load(file=file)
  #genelists
  slist1 <- names(genelists)
  matnorand <- mat[!is.element(mat$source,"Random"),]
  slist2 <- matnorand$signature
  in2not1 <- slist2[!is.element(slist2,slist1)]
  in1not2 <- slist1[!is.element(slist1,slist2)]

  stop <- F
  if(length(in2not1)>0) {
    cat("=========================================================\n")
    cat("In the catalog not in the database\n")
     cat("=========================================================\n")
    print(in2not1)
    stop <- T
  }
  if(stop) browser()
  #if(length(in1not2)>0) {
  #  cat("=========================================================\n")
  #  cat("In the database not in the catalog\n")
  #  cat("=========================================================\n")
  #  #print(in1not2)
  #  stop <- T
  #}
  mat <- mat[mat[,sigset]==1,]
  alist <- mat[is.element(mat$type,"nondirectional"),]
  blist <- mat[is.element(mat$type,"directional"),]

  blist <- blist[order(blist$parent),]
  counter <- 0
  for(parent in blist$parent) {
    temp <- blist[is.element(blist$parent,parent),]
    if(nrow(temp)!=2) {
      cat("missing pair:",parent,"\n")
      counter <- counter+1
    }
  }
  if(counter>0 || stop) {
    browser()
  }
  else return(mat)
}

