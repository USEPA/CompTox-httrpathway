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

  names(mat) <- c("signature","source","subsource","description","gene.list","ngene")
  mat$direction <- "nondirectional"
  mat$type <- "nondirectional"
  mat$parent <- mat$signature
  x <- mat$signature
  for(i in 1:length(x)) {
    sig <- x[i]
    len <- nchar(sig)
    #if(contains(sig,"_UP")) browser()
    if(substr(sig,len-2,len)=="_UP") {
      mat[i,"direction"] <- "up"
      mat[i,"type"] <- "directional"
      mat[i,"parent"] <- substr(sig,1,len-3)
    }
    if(substr(sig,len-2,len)=="_DN") {
      mat[i,"direction"] <- "dn"
      mat[i,"type"] <- "directional"
      mat[i,"parent"] <- substr(sig,1,len-3)
    }
  }
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
  mat[is.element(mat$subsource,class.set),"type"] <- "unidirectional"
  mat[is.element(mat$subsource,class.set),"direction"] <- "both"


  mat.uni <- mat[is.element(mat$type,"unidirectional"),]
  mat.bid <- mat[is.element(mat$type,"bidirectional"),]
  mat.bid <- mat.bid[order(mat.bid$parent),]
  mat <- rbind(mat.bid,mat.uni)
  mat <- mat[order(mat$signature),]

  MsigDB_signatures <- mat
  file <- "../input/signatures/MsigDB_signatures.RData"
  save(MsigDB_signatures,file=file)
}
