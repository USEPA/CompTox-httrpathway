#--------------------------------------------------------------------------------------
#' Summarize the signature catalog
#'
#'--------------------------------------------------------------------------------------
signatureCatalogSummary <- function(sigcatalog="signatureDB_master_catalog 2020-08-31",
                           sigset="screen_large") {
  printCurrentFunction()

  file <- paste0("../input/signatures/",sigcatalog,".xlsx")
  catalog <- read.xlsx(file)
  catalog <- catalog[catalog[,sigset]==1,]

  st.list <- sort(unique(catalog$super_target))
  name.list <- c("super_target","n_signature","n_parent","new_target")
  nst <- length(st.list)
  res <- as.data.frame(matrix(nrow=nst,ncol=length(name.list)))
  names(res) <- name.list
  for(i in 1:nst) {
    st <- st.list[i]
    temp <- catalog[is.element(catalog$super_target,st),]
    n1 <- length(unique(temp$signature))
    n2 <- length(unique(temp$parent))
    res[i,"super_target"] <- st
    res[i,"new_target"] <- st
    res[i,"n_signature"] <- n1
    res[i,"n_parent"] <- n2
  }
  file <- paste0("../input/signatures/",sigcatalog,"_summary.xlsx")
  write.xlsx(res,file)
}
