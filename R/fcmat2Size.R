#--------------------------------------------------------------------------------------
#' Calculate the size of the FCMAT2 files
#' DMEM_6hr_pilot_normal_pe_1
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#' u2os_toxcast_pfas_pe1_normal
#--------------------------------------------------------------------------------------
fcmat2Size <- function() {
  printCurrentFunction()

  ds.list <- c(
    "DMEM_6hr_pilot_normal_pe_1",
    "heparg2d_toxcast_pfas_pe1_normal",
    "mcf7_ph1_pe1_normal_good_pg",
    "u2os_toxcast_pfas_pe1_normal"
  )
  for(dataset in ds.list) {
    file = paste0("../input/fcdata/FCMAT2_",dataset,".RData")
    #print(file)
    load(file=file)
    cat(dataset,"\n")
    cat("   genes: ",ncol(FCMAT2),"\n")
    cat("   samples: ",nrow(FCMAT2),"\n")
    file = paste0("../input/fcdata/CHEM_DICT_",dataset,".RData")
    #print(file)
    load(file=file)
    cat("   chemicals: ",length(unique(CHEM_DICT[,"name"])),"\n")
    browser()

  }
}
