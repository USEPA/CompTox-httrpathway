#' Examine the expression levels of the Corton signature
#'
CortonSignatureExpression = function(do.rawprep=F) {

  printCurrentFunction()

  file = paste0("../input/signatures/signatureDB_genelists.RData")
  print(file)
  load(file=file)

  ryan = sigdb[is.element(sigdb$source,"Ryan"),]
  gl1 = ryan[1,"gene.list"]
  gl2 = ryan[2,"gene.list"]
  gl1g = str_split(gl1,"\\|")[[1]]
  gl2g = str_split(gl2,"\\|")[[1]]
  gl = c(gl1g,gl2g)

  welltype.list <- c("untreated","vehicle_control")
  cutoff=0.95
  dirin = "../input/httr_mcf7_screen/raw_ctrl/"
  dirout=dirin

  if(do.rawprep) {
    for(welltype in welltype.list) {
      file <- paste0(dirin,"httr-ph-i-raw-pl-",welltype,"-v1.tsv")
      print(file)
      mat <- read.delim(file,stringsAsFactors=F)
      x <- unite(mat[,c("plate_id","block_id","pg_id","well_id")],"rowkey")
      mat$rowkey <- x[,1]
      CTRLMAT <- mat
      pmat <- unique(CTRLMAT[,c("plate_id","block_id","pg_id","well_id")])
      x <- unite(pmat,"rowkey")
      rownames(pmat) <- x[,1]
      PLATEMAT <- pmat

      mat <- CTRLMAT[,c("rowkey","probe_id","probe_count")]
      print(dim(mat))
      cat("    start casting\n")
      countmat <- dcast(mat,rowkey~probe_id)
      cat("   casting is done\n")
      rownames(countmat) <- countmat[,"rowkey"]
      countmat <- countmat[,2:ncol(countmat)]
      countmat <- as.matrix(countmat)
      imat <- countmat
      imat[!is.na(imat)] <- 1
      imat[is.na(imat)] <- 0
      cs <- colSums(imat)
      cs <- cs/nrow(countmat)
      countmat <- countmat[,cs>cutoff]
      countmat[is.na(countmat)] <- 0.5
      for(i in 1:nrow(countmat)) countmat[i,] <- 1000000*countmat[i,]/sum(countmat[i,])
      countmat <- log2(countmat)
      COUNTMAT <- countmat
      file <- paste0(dirout,"data_",welltype,"_",cutoff,".RData")
      save(CTRLMAT,PLATEMAT,COUNTMAT,file=file)
    }
  }

  res = as.data.frame(matrix(nrow=length(gl),ncol=3))
  names(res) = c("gene",welltype.list)
  rownames(res) = gl
  res$gene = gl
  for(welltype in welltype.list) {
    cat(welltype,"\n")
    file <- paste0(dirout,"data_",welltype,"_",cutoff,".RData")
    load(file=file)
    x = colMeans(COUNTMAT)
    name.list = names(x)
    y = unlist(str_split(name.list,"_"))
    indx = seq(from=1,to=length(y),by=2)
    z = y[indx]
    names(x) = z
    res[,welltype] = -5
    for(gene in gl) {
      if(is.element(gene,names(x))) {
        res[gene,welltype] = x[gene]
      }
    }
  }
  file = paste0(dirin,"CortonSignature.xlsx")
  write.xlsx(res,file)
  browser()

}

