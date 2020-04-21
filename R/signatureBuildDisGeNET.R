#--------------------------------------------------------------------------------------
#' Build the standard input file for the Ryan signatures
#'
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
signatureBuildDisGeNET <- function(){
  printCurrentFunction()

  file <- "../input/signatures/DisGeNET/gene_disease_raw.txt"
  temp <- readLines(file)
  cat("File read in\n")
  geneid.list <- NULL
  disease.list <- NULL
  for(i in 1:length(temp)) {
    x <- temp[i]
    a <- str_locate_all(pattern="\\[",x)[[1]]
    b <- str_locate_all(pattern="\\]",x)[[1]]
    #cat(i,nrow(temp),"[",nrow(a),"] [",nrow(b),"]\n")
    if(nrow(a)>0 && nrow(b)>0) {
      start <- a[1,1]
      stop <- b[2,1]
      y <- substr(x,start,stop)
      #cat(y,"\n")

      a <- str_locate_all(pattern="\\[",y)[[1]]
      b <- str_locate_all(pattern="\\]",y)[[1]]
      blockg <- substr(y,a[1,1],b[1,1])
      #cat(blockg,"\n")
      geneid <- str_replace_all(blockg,"\\]","")
      geneid <- str_replace_all(geneid,"\\[ncbigene:","")
      geneid <- as.numeric(geneid)
      #cat(geneid,"\n")

      blockd <- substr(y,b[1,1],a[2,1])
      #cat(blockd,"\n")
      disease <- str_replace_all(blockd,"\\] -","")
      disease <- str_replace_all(disease,"\\[","")
      disease <- str_trim(disease)
      #cat(disease,"\n")
      geneid.list <- c(geneid.list,as.numeric(geneid))
      disease.list <- c(disease.list,disease)
    }
    if(i%%1000==0) cat("finished processessing ",i," out of ",length(temp),"\n")
  }
  mat <- as.data.frame(cbind(geneid.list,disease.list,stringsAsFactors=F))
  cat("file processed\n")

  if(!exists("HUMAN.GENES")) {
    file <- "../input/signatures/DisGeNET/gene_info_human.txt"
    gi <- read.table(file,header=F,stringsAsFactors=F,sep="\t",comment="",quote="")
    gi <- gi[gi[,1]==9606,]
    gi <- gi[,c(2:3)]
    names(gi) <- c("geneid","symbol")
    HUMAN.GENES <<- gi
    file <- "../input/signatures/DisGeNET/human_genes.xlsx"
    write.xlsx(HUMAN.GENES,file)
  }
  names(mat) <- c("geneid","disease")
  mat$symbol <- NA
  geneid.list <- unique(mat$geneid)
  for(geneid in geneid.list) {
    if(is.element(geneid,HUMAN.GENES$geneid)) symbol <- HUMAN.GENES[HUMAN.GENES$geneid==geneid,"symbol"]
    else symbol <- "NOGENE"
    mat[mat$geneid==geneid,"symbol"] <- symbol
  }
  mat <- mat[!is.element(mat$symbol,"NOGENE"),]
  mat0 <- mat
  cat("disease gene map created\n")
  disease.list <- unique(mat$disease)
  ndisease <- length(disease.list)

  name.list <- c("signature","set","class","description","gene.list","ngene")
  mat <- as.data.frame(matrix(nrow=ndisease,ncol=length(name.list)))
  names(mat) <- name.list
  mat$signature <- as.character(disease.list)
  mat$set <- "DisGeNET"
  mat$class <- "-"
  mat$description <- as.character(disease.list)
  rownames(mat) <- mat$signature
  for(i in 1:ndisease) {
    disease <- mat[i,"signature"]
    gene.list <- mat0[is.element(mat0$disease,disease),"symbol"]
    mat[disease,"ngene"] <- length(gene.list)
    mat[disease,"gene.list"] <- paste(gene.list,collapse="|")
  }
  cat("final data frame created\n")
  DisGeNET_signatures <- mat
  file <- "../input/signatures/DisGeNET_signatures.RData"
  save(DisGeNET_signatures,file=file)
}
