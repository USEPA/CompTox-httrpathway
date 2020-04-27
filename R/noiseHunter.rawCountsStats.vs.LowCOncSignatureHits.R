#--------------------------------------------------------------------------------------
#' Look for trends in the raw data vs. the low concentraiton signature hits
#'
#' @param dir The directory from which to read all of the raw filesatalog file
#--------------------------------------------------------------------------------------
noiseHunter.rawCountsStats.vs.LowCOncSignatureHits <- function(to.file=F) {
  printCurrentFunction()

  file <- "../output/noiseHunter/noiseHunter.signatureHits.vs.rawCountsStats.xlsx"
  mat <- read.xslx(file)
  if(to.file) {
    fname <- paste0("../output/noiseHunter/noiseHunter raw stats vs low conc hits.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,4),mar=c(4,3,3,2))

  npg <- 48
  if(do.read) {
    dir <- "../output/signature_conc_resp_summary/"
    file <- paste0(dir,"SIGNATURE_CR_screen_large_DMEM_6hr_screen_normal_pe_1_mygsea_0.05_conthits.RData")
    print(file)
    load(file)
    SCR <<- SIGNATURE_CR
    file <- "../output/noiseHunter/httr_mcf7_ph1_metadata.xlsx"
    smap <- read.xlsx(file)
    SMAP <<- smap
  }
  scr <- SCR
  smap <- SMAP
  rownames(smap) <- smap$sample_id
  smap <- smap[is.element(smap$stype,"test sample"),]
  name.list <- c("dtxsid","casrn","name","sample_id","block_id","pg_id",
                 "nhit","nhit.10","nhit.1","nhit.0.1","n_reads_mapd","mapd_frac","bad_probe_count",
                 "n_cov5","n_sig80","top10_prop","gini_coef")
  result <- NULL
  for(pg in 1:npg) {
    cat(pg,"\n")
    smap.pg <- smap[smap$pg_id==pg,]
    smap.pge <- smap.pg[,c("sample_id","well_id","trt_name","stype",
                           "chem_id","dose_level","block_id","pg_id","culture_id",
                           "n_reads","n_reads_mapd","mapd_frac","bad_probe_count",
                           "n_cov5","n_sig80","top10_prop","gini_coef")]
    cid.list <- unique(smap.pg$chem_id)
    res <- as.data.frame(matrix(nrow=length(cid.list),ncol=length(name.list)))
    names(res) <- name.list
    res$sample_id <- cid.list
    res$pg_id <- pg
    res$block_id <- smap.pg[1,"block_id"]
    for(k in 1:length(cid.list)) {
      cid <- cid.list[k]
      temp <- scr[is.element(scr$sample_id,cid),]
      res[k,"casrn"] <- temp[1,"casrn"]
      res[k,"dtxsid"] <- temp[1,"dtxsid"]
      res[k,"name"] <- temp[1,"name"]
      res[k,"sample_id"] <- cid

      temp <- temp[temp$hitcall>0.5,]
      res[k,"nhit"] <- nrow(temp)
      temp <- temp[temp$bmd<10,]
      res[k,"nhit.10"] <- nrow(temp)
      temp <- temp[temp$bmd<1,]
      res[k,"nhit.1"] <- nrow(temp)
      temp <- temp[temp$bmd<0.1,]
      res[k,"nhit.0.1"] <- nrow(temp)

      stemp <- smap.pg[is.element(smap.pg$chem_id,cid),]
      col.list <- names(res)[12:ncol(res)]
      for(col in col.list) res[k,col] <- mean(stemp[,col])
    }
    result <- rbind(result,res)
    #browser()
  }
   write.xlsx(result,file)
}
