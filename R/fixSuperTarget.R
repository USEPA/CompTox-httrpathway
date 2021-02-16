#--------------------------------------------------------------------------------------
#' Replace the super_target values in the signature output file
#' with ones from a new catalog
#'
#' @param dataset The L2fc matrix data set
#' @param dir The directory where the data file lives
#' @param do.read If TRUE, read in FCMAT2 to a gloabal
#' @return
#'
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#' heparg2d_toxcast_pfas_pe1_normal
#' u2os_pilot_pe1_normal_null_pilot_lowconc
#'
#--------------------------------------------------------------------------------------
fixSuperTarget <- function(do.read=T,
                           dataset="u2os_pilot_pe1_normal_null_pilot_lowconc",
                           sigcatalog="signatureDB_master_catalog 2021-02-10",
                           sigset="screen_large",
                           method="fc") {
  printCurrentFunction(paste0(dataset,sigset,method))

  if(do.read) {
    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    MAT <<-SIGNATURE_CR
    file <- paste0("../output/signature_conc_resp_summary/replaced/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    save(SIGNATURE_CR,file=file)
  }
  file <- paste0("../input/signatures/",sigcatalog,".xlsx")
  catalog <- read.xlsx(file)
  catalog <- catalog[catalog[,sigset]==1,]

  mat <- MAT
  parent.list <- unique(catalog$parent)
  cat("mat:",nrow(mat),"\n")
  cat("remove the deprecated signatures\n")
  mat <- mat[is.element(mat$signature,parent.list),]
  cat("mat:",nrow(mat),"\n")
  cat("replace the super_targets\n")
  sig.list <- sort(unique(mat$signature))
  temp <- catalog[is.element(catalog$parent,sig.list),]
  st.list <- sort(unique(temp$super_target))
  for(i in 1:length(st.list)) {
    st <- st.list[i]
    sig.list <- catalog[is.element(catalog$super_target,st),"parent"]
    mat[is.element(mat$signature,sig.list),"super_target"] <- st
    cat(i,length(st.list),":",st,":",length(sig.list),"\n")
  }

  SIGNATURE_CR <- mat
  file <- paste0("../output/signature_conc_resp_summary/fixed/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
  cat("saving ...\n",file,"\n")
  save(SIGNATURE_CR,file=file)
}
