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
signatureScoreMerge <- function(sigset="smallset",
                              sigcatalog="signatureDB_master_catalog 2020-01-31",
                              dataset="DMEM_6hr_pilot_normal_pe_1_RAND10",
                              method="mygsea",
                              nullset="DMEM_6hr_pilot_normal_pe_1_RAND10") {

  printCurrentFunction(paste(dataset,sigset,method))
  starttime = proc.time()

  # get signaturescoremat
  file <- paste0("../output/signature_score_summary/signaturescoremat_",sigset,"_",dataset,"_",method,"_bidirectional.RData")
  print(file)
  load(file)

  sig.list <- unique(signaturescoremat$signature)
  file = paste0("../input/signatures/",sigcatalog,".xlsx")
  annotations <- read.xlsx(file)
  annotations <- annotations[is.element(annotations$signature,sig.list),]
  rownames(annotations) <- annotations$signature
  sig.nondir <- annotations[is.element(annotations$type,"nondirectional"),"signature"]
  sig.updn <- annotations[is.element(annotations$type,"bidirectional"),"signature"]
  par.updn <- unique(annotations[is.element(annotations$type,"bidirectional"),"parent"])

  seta <- signaturescoremat[is.element(signaturescoremat$signature,sig.nondir),]
  setb <- signaturescoremat[is.element(signaturescoremat$signature,sig.updn),]
  setb2 <- NULL
  for(parent in par.updn) {
    sigpar <- annotations[is.element(annotations$parent,parent),]
    sigup <- sigpar[is.element(sigpar$direction,"up"),"signature"]
    sigdn <- sigpar[is.element(sigpar$direction,"dn"),"signature"]

    for(dtxsid in unique(signaturescoremat$dtxsid)) {
      temp1 <- signaturescoremat[is.element(signaturescoremat$dtxsid,dtxsid),]
      temp1up <- temp1[is.element(temp1$signature,sigup),]
      temp1dn <- temp1[is.element(temp1$signature,sigdn),]
      for(conc in unique(temp1up$conc)) {
        temp1upc <- temp1up[temp1up$conc==conc,]
        temp1dnc <- temp1dn[temp1dn$conc==conc,]

        temp <- temp1upc
        temp[1,"signature"] <- parent
        temp[1,"signature_score"] <- temp1upc[1,"signature_score"]-temp1dnc[1,"signature_score"]
        setb2 <- rbind(setb2,temp)
      }
    }
  }
  signaturescoremat <- rbind(seta,setb2)
  file <- paste0("../output/signature_score_summary/signaturescoremat_",sigset,"_",dataset,"_",method,".RData")
  cat("   ",file,"\n")
  save(signaturescoremat,file=file)
}

