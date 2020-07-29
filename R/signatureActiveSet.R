#--------------------------------------------------------------------------------------
#'
#' Create the active set of signatures to make for easier loading
#' @param to.file If TRUE, write plots to a file
#--------------------------------------------------------------------------------------
signatureActiveSet <- function(dataset="DMEM_6hr_screen_normal_pe_1",
                               sigset="screen_large",
                               method = "mygsea",
                               hitcall.threshold=0.5) {
  printCurrentFunction()

  file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
  print(file)
  load(file=file)
  x <- SIGNATURE_CR
  x <- x[x$hitcall>hitcall.threshold,]
  SIGNATURE_CR_ACTIVE <<-x
  file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits active.RData")
  save(SIGNATURE_CR_ACTIVE,file=file)
}

