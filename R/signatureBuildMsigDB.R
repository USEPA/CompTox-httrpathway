#--------------------------------------------------------------------------------------
#' Build the standard input file for the MSigDB signatures
#'
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
signatureBuildMsigDB <- function(){
  printCurrentFunction()

  file <- "../input/signatures/msigdb_PATHWAYS.RData"
  load(file)

  mat <- msigdb_PATHWAYS
  for(i in 1:nrow(mat)) {
    mat[i,"pathway"] <- str_replace_all(mat[i,"pathway"],"_UP","_up")
    mat[i,"pathway"] <- str_replace_all(mat[i,"pathway"],"_DN","_dn")
    mat[i,"pathway"] <- str_replace_all(mat[i,"pathway"],"_DOWN","_dn")
  }

  names(mat) <- c("signature","source","subsource","description","gene.list","ngene")
  mat$direction <- "nondirectional"
  mat$type <- "nondirectional"
  mat$parent <- mat$signature
  x <- mat$signature
  for(i in 1:length(x)) {
    sig <- x[i]
    len <- nchar(sig)
    #if(contains(sig,"_UP")) browser()
    if(substr(sig,len-2,len)=="_up") {
      mat[i,"direction"] <- "up"
      mat[i,"type"] <- "directional"
      mat[i,"parent"] <- substr(sig,1,len-3)
    }
    else if(substr(sig,len-2,len)=="_dn") {
      mat[i,"direction"] <- "dn"
      mat[i,"type"] <- "directional"
      mat[i,"parent"] <- substr(sig,1,len-3)
    }
  }
  cat("==========================================================\n")
  cat("Round 1\n")
  cat("==========================================================\n")
  counter <- 0
  for(i in 1:nrow(mat)) {
    parent <- mat[i,"parent"]
    if(mat[i,"type"]=="directional") {
      temp <- mat[is.element(mat$parent,parent),]
      if(nrow(temp)!=2) {
        counter <- counter+1
        cat("missing pair for parent",parent,":",counter,"\n")
        #print(temp[,1:3])
      }
      mat[i,"type"] <- "nondirectional"
      mat[i,"direction"] <- "nondirectional"
    }
  }
  cat("==========================================================\n")
  cat("Round 2\n")
  cat("==========================================================\n")
  counter <- 0
  for(i in 1:nrow(mat)) {
    parent <- mat[i,"parent"]
    if(mat[i,"type"]=="directional") {
      temp <- mat[is.element(mat$parent,parent),]
      if(nrow(temp)!=2) {
        counter <- counter+1
        cat("missing pair for parent",parent,":",counter,"\n")
        #print(temp[,1:3])
      }
      mat[i,"type"] <- "nondirectional"
      mat[i,"direction"] <- "nondirectional"
    }
  }
  temp <- mat[is.element(mat$type,"directional"),]
  nb <- length(unique(temp$parent))
  cat("directional",nb,"\n")
  name.list <- c("signature","parent","source","type","direction","description","subsource","ngene","gene.list")
  mat <- mat[,name.list]
  class.set <- c("C4: cancer gene neighborhoods","C2: KEGG gene sets",
                 "C2: biocarta gene sets","C2: REACTOME gene sets",
                 "C2: chemical and genetic perturbations","C6: oncogenic gene sets","C7: immunologic signatures",
                 "H: hallmark gene sets","C5: GO biological process","C5: GO cellular component",
                 "c5: GO molecular function","ARCHIVED")
  class.set <- c("C3: transcription factor targets","C4: cancer gene neighborhoods","C2: KEGG gene sets",
                "C2: biocarta gene sets","C2: REACTOME gene sets",
                "C6: oncogenic gene sets","C2: canonical pathways",
                "C5: GO biological process","C5: GO cellular component",
                "c5: GO molecular function","ARCHIVED")
  mat[is.element(mat$subsource,class.set),"type"] <- "nondirectional"
  mat[is.element(mat$subsource,class.set),"direction"] <- "nondirectional"

  mat.non <- mat[is.element(mat$type,"nondirectional"),]
  mat.bid <- mat[is.element(mat$type,"directional"),]
  mat.bid <- mat.bid[order(mat.bid$parent),]
  mat <- rbind(mat.bid,mat.non)
  mat <- mat[order(mat$signature),]
  MsigDB_signatures <- mat
  file <- "../input/signatures/MsigDB_signatures.RData"
  save(MsigDB_signatures,file=file)
}
