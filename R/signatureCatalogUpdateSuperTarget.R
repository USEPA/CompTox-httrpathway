#--------------------------------------------------------------------------------------
#' Update the super_target values in the signature catalog
#'
#'--------------------------------------------------------------------------------------
signatureCatalogUpdateSuperTarget <- function(sigcatalog.in="signatureDB_master_catalog 2020-10-29",
                                              sigcatalog.out="signatureDB_master_catalog 2020-08-31",
                                              mapfile="signatureDB_master_catalog summary in",
                                              sigset="screen_large") {
  printCurrentFunction()

  file <- paste0("../input/signatures/",sigcatalog.in,".xlsx")
  catalog.in <- read.xlsx(file)
  catalog.out = catalog.in
  file <- paste0("../input/signatures/",mapfile,".xlsx")
  smap <- read.xlsx(file)

  for(i in 1:nrow(smap)) {
    stin = smap[i,"super_target"]
    stout = smap[i,"new_target"]
    if(stout!=stin) {
      cat(stin,":",stout,"\n")
      catalog.out[is.element(catalog.out$super_target,stin),"super_target"] = stout
    }
  }
  file <- paste0("../input/signatures/",sigcatalog.out,".xlsx")
  write.xlsx(catalog.out,file)
  signatureCatalogSummary(sigcatalog=sigcatalog.out,sigset=sigset)
}
