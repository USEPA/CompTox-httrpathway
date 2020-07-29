#' understand why the signatures cores differ by number of signatures being calculated
#'
#' Retrieves signature cutoffs for a given null dataset.
signatureScoreDebug = function(dataset="DMEM_6hr_pilot_normal_pe_1",
                               nullset="DMEM_6hr_pilot_normal_pe_1_RAND1000",
                               method="mygsea",
                               signature="CMAP fulvestrant 1e-06 100 8242 100",
                               name="Fulvestrant"){
  printCurrentFunction(paste(dataset,":",nullset))

  sigset.list <- c("pilot_tiny","test20A","test20AB")
  cat("-------------------------------------------------\n")
  cat("nullset\n")
  cat("-------------------------------------------------\n")
  for(sigset in sigset.list) {
    file <- paste0("../output/signature_score_summary/signaturescoremat_",sigset,"_",nullset,"_",method,".RData")
    load(file)
    temp <- signaturescoremat[is.element(signaturescoremat$signature,signature),]
    mval <- mean(temp$signature_score)
    cat(sigset,nrow(temp),mean(temp$signature_score),sd(temp$signature_score),"\n")
  }

  cat("-------------------------------------------------\n")
  cat("dataset\n")
  cat("-------------------------------------------------\n")
  for(sigset in sigset.list) {
    file <- paste0("../output/signature_score_summary/signaturescoremat_",sigset,"_",dataset,"_",method,".RData")
    load(file)
    temp <- signaturescoremat[is.element(signaturescoremat$signature,signature),]
    mval <- mean(temp$signature_score)
    cat(sigset,nrow(temp),mean(temp$signature_score),sd(temp$signature_score),"\n")
  }

}
