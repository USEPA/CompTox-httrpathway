#' Merge the up and down halves of the directional pathway data
#'
#' @param sigset Name of the signature set.
#' @param sigcatalog Name of the catalog file
#' @param method Pathway scoring method in c("fc", "gsva", "gsea")
#' @param signaturescoremat dataframe returned by the upstream signatureScore function
#' @param sigdbgenelist full path to signature DB gene list file; default is repo version
#'
#' @importFrom stringi stri_replace
#' @return signaturescoremat dataframe
#' @export signatureScoreMerge
signatureScoreMerge <- function(sigset,
                                sigcatalog="../inst/extdata/signatureDB_master_catalog_2022-05-16.xlsx",
                                method,
                                signaturescoremat,
                                sigdbgenelist="../inst/extdata/signatureDB_genelists.RDS") {


  printCurrentFunction(paste(sigset,method))
  starttime = proc.time()

  sig.list <- unique(signaturescoremat$signature)

  annotations <- signatureCatalogLoader(sigset, sigcatalog, sigdbgenelist)
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

    # [LJE 12/22/23]
    # Adding code to drop the unpaired signatures here
    # The better place to do this would be in signatureScore.R but it's not clear the full signature catalog is actually loaded there?
    sig_parents <- unique(setb[, c("signature", "parent")])
    stopifnot(sum(duplicated(sig_parents$signature)) == 0)
    cat("[Debug:]", nrow(sig_parents), "directional signatures were scored corresponding to", length(unique(sig_parents$parent)), "parents.\n")
    # Check if there are any parents with only one half here:
    n_sig_parents <- table(sig_parents$parent)
    stopifnot(all(n_sig_parents %in% c(1,2)))
    paired_sigs <- names(n_sig_parents)[n_sig_parents == 2]
    unpaired_sigs <- names(n_sig_parents)[n_sig_parents == 1]
    if(length(unpaired_sigs) > 0) {
      cat("Dropping", length(unpaired_sigs), "bidirectional signatures missing up or dn half, e.g.", paste(head(unpaired_sigs), collapse = ", "), "\n")
      setb <- setb[setb$parent %in% paired_sigs, ]
      # Make sure what's left are the paired parents
      stopifnot(all(paired_sigs %in% setb$parent))
      stopifnot(all(setb$parent %in% paired_sigs))
      # Now the assumption that all bidirectional signatures have an up and a dn half is valid
      cat(length(unique(setb$signature)), "bidirectional signatures",
          "corresponding to", length(unique(setb$parent)), "parent signatures remain.\n")
    }
    # Merge the type, parent, and direction columns from annotations into setb...
    stopifnot(all(setb$signature %in% annotations$signature))
    setb_annot <- merge(setb, annotations[, c("signature", "type", "direction", "parent")],
                        by = "signature", all.x = T, all.y = F, sort = F, suffixes = c("", "_orig"))
    stopifnot(nrow(setb) == nrow(setb_annot))
    # NOTE: Setting sort = F in merge doesn't mean the rows come out in the same order...
    stopifnot(all(sort(setb$index) == sort(setb_annot$index)))
    stopifnot(all(setb_annot$parent == setb_annot$parent_orig))
    stopifnot(all(setb_annot$type == "directional"))
    stopifnot(all(setb_annot$direciton %in% c("dn", "up")))
    setb <- setb_annot
    # Resume original code...

    cat("> order the table\n")
    setc <- setb[order(setb$index),]

  cat("> split the table\n")
  # [LJE 12/22/23]
  # The problem was here - it assumes there's still a matching number of up and dn signatures, but that's not the case
  # Even when this doesn't throw an error, if there was an imbalance the up and dn signatures may not match up...
  # Modified the code below to split on the actual direction column instead of assuming the sort will guarantee this
  # top <- nrow(setc)
  #  index1 <- seq(from=1,to=(top-1),by=2)
  #  index2 <- index1+1

  #  setc1 <- setc[index1,] #dn
  #  setc2 <- setc[index2,] #up
  setc1 <- setc[setc$direction == "dn", ]
  setc2 <- setc[setc$direction == "up", ]
    # Adding checks here to make sure setc1 is actually all dn signatures and setc2 is actually all up signatures
    if(!all(grepl("dn$", setc1$signature))) {
      warning("setc1 should be all signatures ending in \"dn\" but contains other signatures.\n")
      stop("setc1 should be all signatures ending in \"dn\" but contains other signatures.\n")
    }
    if(!all(grepl("up$", setc2$signature))) {
      warning("setc2 should be all signatures ending in \"up\" but contains other signatures.\n")
      stop("setc2 should be all signatures ending in \"up\" but contains other signatures.\n")
    }
    # Resume original code...
    cat("> rows in setc1 and set c2 before double check",nrow(setc1),nrow(setc2),"\n")
    setc1$index <- paste(setc1$dtxsid,setc1$sample_id,setc1$conc,setc1$parent)
    setc2$index <- paste(setc2$dtxsid,setc2$sample_id,setc2$conc,setc2$parent)
    mask1 <- duplicated(setc1$index)
    mask2 <- duplicated(setc2$index)
    cat("duplicates in Set 1:",sum(mask1),length(mask1),"\n")
    cat("duplicates in Set 2:",sum(mask2),length(mask2),"\n")

    # [LJE 12/22/23]
    # A few final checks here - there should be no duplicates and the indexes should match completely between setc1 and setc2
    if(sum(mask1) > 0 || sum(mask2) > 0) {
      warning("There are still duplicate index values in setc1 or setc2.\n")
      stop("There are still duplicate index values in setc1 or setc2.\n")
    }
    if(nrow(setc1) == nrow(setc2)) {
      cat("Split bidirectional signature score matrix into dn and up halves with", nrow(setc1), "rows each.\n")
    } else {
      warning("setc1 and setc2 have different numbers of rows.\n")
      stop("setc1 and setc2 have different numbers of rows.\n")
    }
    # Make sure the indexes match up now that they're named by parent
    if(!all(sort(setc1$index) == sort(setc2$index))) {
      warning("There is still mismatching index values in setc1 vs setc2")
      stop("There is still mismatching index values in setc1 vs setc2")
    }

    rownames(setc1) <- setc1$index
    rownames(setc2) <- setc2$index

    index <- setc1$index
    index <- index[is.element(index,setc2$index)]
    setd1 <- setc1[is.element(setc1$index,index),]
    setd2 <- setc2[is.element(setc2$index,index),]
    setd1 <- setd1[index,]
    setd2 <- setd2[index,]

    cat("> rows in setc1 and set c2 after double check",nrow(setd1),nrow(setd2),"\n")
    # [LJE 12/22/23] One last sanity check...
    stopifnot(all(setd1$index == setd2$index))

    pc1 <- unique(setc1$parent)
    pd1 <- unique(setd1$parent)
    plost <- pc1[!is.element(pc1,pd1)]
    if(length(plost)>0) {
      warning("> unmatched up-down signatures\n")
      for(i in 1:length(plost)) warning(plost[i],"\n")
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
  return(signaturescoremat)
}

