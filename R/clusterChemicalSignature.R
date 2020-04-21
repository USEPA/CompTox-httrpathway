#--------------------------------------------------------------------------------------
#' Run the clustering of chemicals and signatures
#'
#' @param min.ngene Signatures will only be saved if the number of genes is >= this value
#' @param max.ngene Signatures will only be saved if the number of genes is <= this value
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
clusterChemicalSignature <- function(do.load=F,
                                     sigset="pilot_large_all_CMAP",
                                     sigcatalog="signatureDB_master_catalog 2020-03-10",
                                     dataset="DMEM_6hr_screen_normal_pe_1",
                                     method="mygsea",
                                     cutoff=1,
                                     k=5){
  printCurrentFunction()

  if(do.load) {
    #file <- paste("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData",sep="")
    #print(file)
    #load(file=file)
    #signaturecr <<- SIGNATURE_CR
    annotations <<- signatureCatalogLoader(sigset,sigcatalog)
  }

  mat <- signaturecr
  #mat <- mat[1:100000,]
  mat[mat$hitcall<0.5,"bmd"] <- 1000
  res <- reshape2::dcast(mat,name~signature,value.var="bmd",fun.aggregate = mean,fill=1000)
  rownames(res) <- res[,1]
  res <- res[,2:ncol(res)]
  res[res>cutoff] <- 1000
  res <- 3-log10(res)
  res[res>0] <- 1
  signatureChemicalBMDmatrix <- res
  file <- paste0("../output/signature_corr/signatureChemicalBMDmatrix ",cutoff," ",k,".RData")
  save(signatureChemicalBMDmatrix,file=file)

  rs <- rowSums(res)
  cs <- colSums(res)

  cat("start clustering\n")
  chem.clusters <- kmeans(res,centers=k)
  sig.clusters <- kmeans(t(res),centers=k)

  x <- as.data.frame(chem.clusters$cluster)
  x <- cbind(x,x)
  x[,1] <- rownames(x)
  names(x) <- c("chemical","cluster")
  rs2 <- rs[rownames(x)]
  x$nhit <- rs2

  file <- "../input/chemicals/screen_chemicals_target_annoations.xlsx"
  temp <- read.xlsx(file)
  rownames(temp) <- temp$name
  temp <- temp[x$chemical,]
  x$target <- temp$target
  chem.cluster.mat <- x
  file <- paste0("../output/signature_corr/chemicalClusters ",cutoff," ",k,".xlsx")
  write.xlsx(chem.cluster.mat,file)

  x <- as.data.frame(sig.clusters$cluster)
  x <- cbind(x,x)
  x[,1] <- rownames(x)
  names(x) <- c("signature","cluster")
  cs2 <- cs[rownames(x)]
  x$nhit <- cs2

  temp <- unique(annotations[,c("parent","super_target")])
  rownames(temp) <- temp$parent
  temp <- temp[x$signature,]
  x$super_target <- temp$super_target

  sig.cluster.mat <- x
  file <- paste0("../output/signature_corr/signatureClusters ",cutoff," ",k,".xlsx")
  write.xlsx(sig.cluster.mat,file)
}
