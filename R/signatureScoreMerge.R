#' Merge the up and down halves of the pathway data
#'
#' @param sigset Name of the signature set.
#' @param sigcatlog Nmae of the catalog file
#' @param dataset Name of the data set.
#' @param method Pathway scoring method in c("fc", "gsva", "mygsea")
#' @param nullset Name of the null data set.
#'
#' @import data.table
#' @import parallel
#' @import openxlsx
#'
#' @return nothing
#' @export
signatureScoreMerge <- function(sigset="screen_large",
                                sigcatalog="signatureDB_master_catalog 2020-04-04",
                                dataset="DMEM_6hr_screen_normal_pe_1_RAND1000",
                                method="mygsea",
                                nullset="DMEM_6hr_screen_normal_pe_1_RAND1000") {

  printCurrentFunction(paste(dataset,sigset,method))
  starttime = proc.time()

  cat("> get signaturescoremat\n")
  file <- paste0("../output/signature_score_summary/signaturescoremat_",sigset,"_",dataset,"_",method,"_directional.RData")
  print(file)
  load(file)

  sig.list <- unique(signaturescoremat$signature)

  annotations <- signatureCatalogLoader(sigset,sigcatalog)
  annotations <- annotations[is.element(annotations$signature,sig.list),]
  rownames(annotations) <- annotations$signature
  sig.nondir <- annotations[is.element(annotations$type,c("nondirectional","unidirectional")),"signature"]
  sig.updn <- annotations[is.element(annotations$type,"directional"),"signature"]
  par.updn <- unique(annotations[is.element(annotations$type,"directional"),"parent"])

  cat("> split the directional from the nondirectional\n")
  seta <- signaturescoremat[is.element(signaturescoremat$signature,sig.nondir),]
  setb <- signaturescoremat[is.element(signaturescoremat$signature,sig.updn),]

  cat(" seta:",dim(seta),"\n")
  cat(" setb:",dim(setb),"\n")
  seta$parent <- seta$signature

  if(nrow(setb)>0) {
    cat("> split the up and down from the signature to get the parent\n")
    x <- setb$signature
    y <- stri_replace(x,replacement="",mode="last",fixed=" up")
    y <- stri_replace(y,replacement="",mode="last",fixed=" dn")
    y <- stri_replace(y,replacement="",mode="last",fixed="_up")
    y <- stri_replace(y,replacement="",mode="last",fixed="_dn")
    setb$parent <- y
    setb$index <- paste(setb$dtxsid,setb$sample_id,setb$conc,setb$signature)
    cat("> order the table\n")
    setc <- setb[order(setb$index),]

  cat("> split the table\n")
  top <- nrow(setc)
    index1 <- seq(from=1,to=(top-1),by=2)
    index2 <- index1+1

    setc1 <- setc[index1,] #dn
    setc2 <- setc[index2,] #up
    cat("> rows in setc1 and set c2 before double check",nrow(setc1),nrow(setc2),"\n")
    setc1$index <- paste(setc1$dtxsid,setc1$sample_id,setc1$conc,setc1$parent)
    setc2$index <- paste(setc2$dtxsid,setc2$sample_id,setc2$conc,setc2$parent)
    mask1 <- duplicated(setc1$index)
    mask2 <- duplicated(setc2$index)
    cat("duplicates in Set 1:",sum(mask1),length(mask1),"\n")
    cat("duplicates in Set 2:",sum(mask2),length(mask2),"\n")
    #browser()
    rownames(setc1) <- setc1$index
    rownames(setc2) <- setc2$index

    index <- setc1$index
    index <- index[is.element(index,setc2$index)]
    setd1 <- setc1[is.element(setc1$index,index),]
    setd2 <- setc2[is.element(setc2$index,index),]
    setd1 <- setd1[index,]
    setd2 <- setd2[index,]

    cat("> rows in setc1 and set c2 after double check",nrow(setd1),nrow(setd2),"\n")

    pc1 <- unique(setc1$parent)
    pd1 <- unique(setd1$parent)
    plost <- pc1[!is.element(pc1,pd1)]
    if(length(plost)>0) {
      cat("> unmatched up-down signatures\n")
      for(i in 1:length(plost)) cat(plost[i],"\n")
      browser()
    }

    sete <- setd2
    sete$signature_score <- setd2$signature_score - setd1$signature_score
    sete$signature <- sete$parent
    sete <- sete[,names(seta)]
    signaturescoremat <- rbind(seta,sete)
  }
  else {
    signaturescoremat <- seta
  }
  file <- paste0("../output/signature_score_summary/signaturescoremat_",sigset,"_",dataset,"_",method,".RData")
  cat("   ",file,"\n")
  save(signaturescoremat,file=file)
}

