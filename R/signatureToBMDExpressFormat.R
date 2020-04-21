library(org.Hs.eg.db)
#' Map the signatures ot the BMDExpress format
#'
#' @param sigset Name of the signature set.
#' @param sigcatlog Nmae of the catalog file
#' @return the trimmed signature table
#' @export
signatureToBMDExpressFormat <- function(sigset="pilot_tiny",
                                   sigcatalog="signatureDB_master_catalog 2020-01-31") {

  printCurrentFunction()

  file = paste0("../input/signatures/",sigcatalog,".xlsx")
  mat <- read.xlsx(file)
  file = paste0("../input/signatures/signatureDB_genelists.RData")
  load(file=file)
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

  mat <- mat[mat[,sigset]==1,]
  x <- org.Hs.egSYMBOL2EG
  mapped_genes <- mappedkeys(x)
  xx <- as.list(x[mapped_genes])
  map <- unlist(xx)
  sigout <- NULL
  counter <- 1
  for(sig in mat$signature) {
    gl <- genelists[sig][[1]]
    gid <- as.numeric(as.vector(map[gl]))
    gid <- gid[!is.na(gid)]
    ngene <- length(gid)
    temp <- as.data.frame(matrix(nrow=ngene,ncol=3))
    names(temp) <- c("pathid","pathname","EntrezGeneID")
    counter <- counter+1

    temp[,1] <- paste0("CCTE_HTTR_",counter)
    temp[,2] <- sig
    temp[,3] <- gid
    sigout <- rbind(sigout,temp)
  }

  file <- paste0("../input/signatures/BMDEXpress pathway input ",sigset,".txt")
  write.table(sigout,file=file,sep="\t",quote=F,row.names=F)
}

