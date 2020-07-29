#--------------------------------------------------------------------------------------
#'
#' Get the signature ranks for chemicals
#--------------------------------------------------------------------------------------
signatureRankScreen <- function(do.load=F,
                                dataset="DMEM_6hr_screen_normal_pe_1",
                                sigset="screen_large",
                                sigcatalog="signatureDB_master_catalog 2020-04-05",
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

  file <- paste0("../input/signatures/signature.target.match ",Sys.Date()," manual.xlsx")
  print(file)
  tmap <- read.xlsx(file)
  tmap <- tmap[tmap$match==1,]

  #dtxsid.list <- unique(smap[smap$block_id==3,"dtxsid"])
  deseq <- deseq[is.element(deseq$dtxsid,dtxsid.list),]
  chems <- chems[is.element(chems$dtxsid,dtxsid.list),]
  nchem <- nrow(chems)

  drank <- abs(deseq$hitcall * log10(deseq$bmd))
  drank <- deseq$hitcall * log10(deseq$bmd)
  deseq$drank <- drank
  deseq <- deseq[order(deseq$drank,decreasing=F),]

  #annotations <- signatureCatalogLoader(sigset,sigcatalog)

  name.list <- c("dtxsid","casrn","name","block_id","pg_id","target",
                 "deseq2.first.signature","deseq2.first.signature.score","deseq2.first.signature.class",
                 "deseq.first.ontarget.signature",
                 "deseq.first.ontarget.signature.class",
                 "deseq.first.ontarget.signature.rank",
                 "deseq.first.ontarget.signature.score",
                 "deseq.first.ontarget.signature.pod")
  result <- as.data.frame(matrix(nrow=nchem,ncol=length(name.list)))
  names(result) <- name.list
  result$dtxsid <- chems$dtxsid
  rownames(result) <- result$dtxsid
  nsig <- length(unique(deseq$signature))
  for(dtxsid in dtxsid.list) {
    temp <- deseq[is.element(deseq$dtxsid,dtxsid),]
    result[dtxsid,"name"] <- temp[1,"name"]
    result[dtxsid,"casrn"] <- temp[1,"casrn"]
    x <- smap[is.element(smap$dtxsid,dtxsid),]
    result[dtxsid,"block_id"] <- x[1,"block_id"]
    result[dtxsid,"pg_id"] <- x[1,"pg_id"]
    target <- chems[dtxsid,"target"]
    result[dtxsid,"target"] <- target
    cat(result[dtxsid,"name"],":",target,"\n")
    target.list <- NULL
    if(target!="-") target.list <- str_split(target,"\\|")[[1]]
    if(!is.null(target.list)) {
      result[dtxsid,"deseq2.first.signature"] <- temp[1,"signature"]
      result[dtxsid,"deseq2.first.signature.score"] <- temp[1,"drank"]
      result[dtxsid,"deseq2.first.signature.class"] <- temp[1,"super_target"]
      mask <- target.match(temp$super_target,target.list,tmap)
      #mask[is.element(temp$super_target,target.list)] <- 1
      if(max(mask)>0) {
        index <- min(which(mask==1))
        result[dtxsid,"deseq.first.ontarget.signature"] <- temp[index,"signature"]
        result[dtxsid,"deseq.first.ontarget.signature.class"] <- temp[index,"super_target"]
        result[dtxsid,"deseq.first.ontarget.signature.rank"] <- index / nsig
        result[dtxsid,"deseq.first.ontarget.signature.score"] <- temp[index,"drank"]
        result[dtxsid,"deseq.first.ontarget.signature.pod"] <- temp[index,"bmd"]
      }
    }
  }
  file <- paste0("../output/signature_rank/signature_rank_",dataset,"_",sigset,"_",method,".xlsx")
  write.xlsx(result,file)
}
target.match <- function(x,target.list,tmap) {
  mask <- vector(length=length(x),mode="integer")
  mask[] <- 0
  for(i in 1:length(x)) {
    targets <- tmap[is.element(tmap$signature,x[i]),"chemical"]
    if(sum(is.element(targets,target.list))>0) {
      mask[i] <- 1
    }
  }
  return(mask)
}
