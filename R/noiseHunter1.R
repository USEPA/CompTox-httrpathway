#--------------------------------------------------------------------------------------
#' Do the first step in finding the noise data - chemicasl that have many signatures
#' with BMD at low concentrations
#'
#' @param do.load If TRUE, load the data into a global
#' @param dataset The name of the data set to analyze
#' @param sigset The signature set to use
#' @param sigcatalog The signature catalog to use
#' @param method The signature scoring method
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
noiseHunter1 <- function(do.load=F,
                         dataset="DMEM_6hr_screen_normal_pe_1",
                         sigset="pilot_large_all_CMAP",
                         sigcatalog="signatureDB_master_catalog 2020-03-10",
                         method="mygsea",
                         cutoff=0.25){
  printCurrentFunction()

  if(do.load) {
    file <- paste("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData",sep="")
    print(file)
    load(file=file)
    signaturecr <<- SIGNATURE_CR
    annotations <<- signatureCatalogLoader(sigset,sigcatalog)

    file <- paste0("../input/chemicals/",dataset,"_chemical_map.xlsx")
    cmap <<- read.xlsx(file)
  }
  cmap.small <- unique(cmap[,c("sample_id","pg_id","block_id")])
  rownames(cmap.small) <- cmap.small$sample_id
  file <- "../input/chemicals/screen_chemicals_target_annoations.xlsx"
  targets <- read.xlsx(file)
  rownames(targets) <- targets$dtxsid
  mat <- signaturecr[signaturecr$hitcall>0.5,]
  mat <- mat[mat$bmd<cutoff,]

  annot <- unique(annotations[,c("parent","super_target")])
  names(annot) <- c("signature","super_target")
  rownames(annot) <- annot$signature
  rownames(annotations) <- annotations$signature
  res <- unique(mat[,c("signature","signature"),])
  res$super_target <- "-"
  res$nsig <- 0
  res$mean_hitcall <- -1
  res$mean_top_over_cutoff <- -1
  res <- res[!is.na(res$signature),]
  for(i in 1:nrow(res)) {
    signature <- res[i,"signature"]
    res[i,"super_target"] <- annot[signature,"super_target"]
    temp <- mat[is.element(mat$signature,signature),]
    res[i,"nsig"] <- nrow(temp)
    res[i,"mean_hitcall"] <- mean(temp$hitcall)
    res[i,"mean_top_over_cutoff"] <- mean(temp$top_over_cutoff)
  }
  res <- res[,2:ncol(res)]
  file <- paste0("../output/noiseHunter/noiseHunter1.signatures ",cutoff,".xlsx")
  write.xlsx(res,file)

  res <- unique(mat[,c("dtxsid","casrn","name","sample_id")])
  res$target <- "-"
  res$pg_id <- -1
  res$block_id <- -1
  res$nsig <- 0
  res$mean_hitcall <- -1
  res$mean_top_over_cutoff <- -1
  res <- res[!is.na(res$dtxsid),]
  for(i in 1:nrow(res)) {
    sample_id <- res[i,"sample_id"]
    dtxsid <- res[i,"dtxsid"]
    res[i,"target"] <- targets[dtxsid,"target"]
    res[i,"pg_id"] <- cmap.small[sample_id,"pg_id"]
    res[i,"block_id"] <- cmap.small[sample_id,"block_id"]
    temp <- mat[is.element(mat$sample_id,sample_id),]
    res[i,"nsig"] <- nrow(temp)
    res[i,"mean_hitcall"] <- mean(temp$hitcall)
    res[i,"mean_top_over_cutoff"] <- mean(temp$top_over_cutoff)
  }
  file <- paste0("../output/noiseHunter/noiseHunter1.chemicals ",cutoff,".xlsx")
  write.xlsx(res,file)
}
