#--------------------------------------------------------------------------------------
#'
#' Get the signature ranks for chemicals
#--------------------------------------------------------------------------------------
signaturePrepTargetMatch <- function(do.load=F,
                                dataset="DMEM_6hr_screen_normal_pe_1",
                                sigset="screen_large",
                                method="mygsea") {
  printCurrentFunction()

  file <- "../input/chemicals/screen_chemicals_target_annoations.xlsx"
  chems <- read.xlsx(file)
  chems <- chems[order(chems$name),]
  dtxsid.list <- chems$dtxsid
  nchem <- length(dtxsid.list)
  rownames(chems) <- chems$dtxsid

  file <- paste0("../input/chemicals/",dataset,"_chemical_map.xlsx")
  print(file)
  smap <- read.xlsx(file)

  if(do.load) {
    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits hits.RData")
    print(file)
    load(file)
    SIGNATURE_CR <<- SIG_CONC_RESPONSE_ACTIVE
  }
  deseq <- SIGNATURE_CR

  sig.targets <- unique(deseq$super_target)
  chem.targets <- unique(chems$target)
  x <- str_split(chem.targets,"\\|")
  y <- sort(unique(unlist(x)))
  result <- expand.grid(signature=sig.targets,chemical=y)
  result$match <- 0

  file <- paste0("../input/signatures/signature.target.match ",Sys.Date(),".xlsx")
  write.xlsx(result,file)
}
