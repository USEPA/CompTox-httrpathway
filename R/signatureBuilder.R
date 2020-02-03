#--------------------------------------------------------------------------------------
#' Create the fiels needed for the signature calculations
#'
#' @param catalog.file This is the name of the catalog file
#' @param pathsetname THis is the name of the pathway set and will be used for the
#' name of the output file
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
signatureBuilder = function(nrandom=1000){
  printCurrentFunction()
  load("../input/signatures/MsigDB_signatures.RData")
  load("../input/signatures/Ryan_signatures.RData")
  load("../input/signatures/Bioplanet_signatures.RData")
  load("../input/signatures/CMAP_signatures.RData")

  sigdb <- rbind(Bioplanet_signatures,MsigDB_signatures,Ryan_signatures,CMAP_signatures)
  rownames(sigdb) <- sigdb$signature

  name.list <- c("signature","parent","source","subsource","type","direction","ngene","description","gene.list")
  sigdb <- sigdb[,name.list]

  x <- sigdb$gene.list
  y <- paste(x,collapse="|")
  genedist <- strsplit(y, "\\|")[[1]]

  nmean <- floor(mean(sigdb$ngene))
  nsd <- floor(sd(sigdb$ngene))
  nsd <- 50

  sigrand <- as.data.frame(matrix(nrow=nrandom,ncol=ncol(sigdb)))
  names(sigrand) <- names(sigdb)
  sigrand$source <- "Random"
  sigrand$type <- "unidirectional"
  sigrand$direction <- "both"
  sigrand$description <- "Random"
  sigrand$subsource <- "-"
  for(i in 1:nrandom) {
    sigrand[i,"signature"] <- paste0("Random_",i)
    sigrand[i,"parent"] <- paste0("Random_",i)
    ngene <- floor(rnorm(1,nmean,nsd))
    if(ngene<10) ngene <- 10
    gene.list <- sample(genedist,ngene)
    sigrand[i,"ngene"] <- ngene
    sigrand[i,"gene.list"] <- paste(gene.list,collapse="|")
  }
  sigdb <- rbind(sigdb,sigrand)

  genelists = strsplit(sigdb$gene.list, "\\|")
  print(length(genelists))
  names(genelists) = sigdb$signature

  catalog <- sigdb[,1:8]
  catalog$target_class <- "-"
  catalog$super_target <- "-"
  catalog$include0 <- 1
  catalog$set1 <- 0
  catalog$set2 <- 0
  catalog$set3 <- 0
  catalog$set4 <- 0
  catalog$set5 <- 0

  file <- "../input/signatures/CMAP refchemdb output.xlsx"
  refchem <- read.xlsx(file)
  rownames(refchem) <- refchem$description
  refchem <- refchem[!is.na(refchem$target),]
  for(i in 1:nrow(refchem)) {
    desc <- refchem[i,"description"]
    target <- refchem[i,"target"]
    catalog[is.element(catalog$description,desc),"target"] <- target
  }

  file <- "../input/signatures/signatureDB.RData"
  save(sigdb,file=file)
  file <- paste0("../input/signatures/signatureDB_genelists.RData")
  save(genelists,file=file)

  file <- paste0("../input/signatures/signatureDB_master_catalog.xlsx")
  write.xlsx(catalog,file)
}
