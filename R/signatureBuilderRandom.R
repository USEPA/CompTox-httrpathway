library(parallel)
#--------------------------------------------------------------------------------------
#' Add the random gene sets to the signature files
#'
#' @param nrandom Number of random gene sets
#' @param mc.cores The number of cores to use in parallel
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
signatureBuilderRandom = function(nrandom=1000,mc.cores=1){
  printCurrentFunction()

  file <- paste0("../input/signatures/signatureDB_genelists no rand.RData")
  load(file=file)
  # genelists

  file <- paste0("../input/signatures/signatureDB_master_catalog no rand.xlsx")
  catalog <- read.xlsx(file)
  file <- "../input/signatures/signatureDB no rand.RData"
  load(file=file)
  #sigdb
  temp <- sigdb[sigdb$ngene<=500,]
  temp <- temp[temp$ngene>=10,]

  cat("> files read in:",nrow(temp)," signatures \n")
  mask <- vector(mode="integer",length=nrow(temp))
  mask[] <- 0
  index <- sample(seq(from=1,to=nrow(temp),by=1),5000)
  mask[index] <- 1
  temp <- temp[mask==1,]
  fl <- temp$gene.list
  nn <- length(fl)
  pairs <- NULL
  pairs <- mclapply(fl,FUN=pairit,mc.cores=mc.cores)

  pairs <- unlist(pairs)
  cat("> pairs created: ",length(pairs),"\n")

  nmean <- floor(mean(temp$ngene))
  nsd <- floor(sd(temp$ngene))
  nsd <- 50

  sigrand <- as.data.frame(matrix(nrow=nrandom,ncol=ncol(temp)))
  names(sigrand) <- names(temp)
  sigrand$source <- "Random"
  sigrand$type <- "nondirectional"
  sigrand$direction <- "both"
  sigrand$description <- "Random"
  sigrand$subsource <- "-"
  for(i in 1:nrandom) {
    sigrand[i,"signature"] <- paste0("Random_",i)
    sigrand[i,"parent"] <- paste0("Random_",i)
    ngene <- floor(rnorm(1,nmean,nsd))
    if(ngene<10) ngene <- 10
    ngene0 <- ngene/2
    pair.list <- sample(pairs,ngene0)
    gene.list <- unlist(strsplit(pair.list,"\\|"))
    gene.list <- unique(gene.list)
    sigrand[i,"ngene"] <- ngene
    sigrand[i,"gene.list"] <- paste(gene.list,collapse="|")
  }
  sigdb <- rbind(sigdb,sigrand)
  cat("> random lists created:",nrow(sigrand),nrow(sigdb),"\n")
  genelists = strsplit(sigdb$gene.list, "\\|")
  cat("> length of gene lists:",length(genelists),"\n")
  names(genelists) = sigdb$signature

  file <- "../input/signatures/signatureDB.RData"
  save(sigdb,file=file)
  file <- paste0("../input/signatures/signatureDB_genelists.RData")
  save(genelists,file=file)
  rcatalog <- sigdb[is.element(sigdb$source,"Random"),1:8]
  rcatalog$target_class <- "Random"
  rcatalog$super_target <- "Random"
  rcatalog$include0 <- 1
  rcatalog$set1 <- 0
  rcatalog$set2 <- 0
  rcatalog$set3 <- 0
  rcatalog$set4 <- 0
  rcatalog$set5 <- 0
  catalog <- rbind(catalog,rcatalog)
  cat("> new catalog created:",nrow(catalog),"\n")
  file <- paste0("../input/signatures/signatureDB_master_catalog.xlsx")
  write.xlsx(catalog,file)
  cat("> file written out\n")
}
pairit <- function(x) {
  gl <- sort(strsplit(x,"\\|")[[1]])
  cm <- t(combn(gl,2))
  cmp <- paste0(cm[,1],"|",cm[,2])
  return(cmp)
}
