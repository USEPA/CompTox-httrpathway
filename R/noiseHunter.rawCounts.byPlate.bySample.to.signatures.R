#--------------------------------------------------------------------------------------
#' Look for trends in the raw data by plate within plate group
#'
#' @param dir The directory from which to read all of the raw filesatalog file
#--------------------------------------------------------------------------------------
noiseHunter.rawCounts.byPlate.bySample.to.signatures <- function(do.read=F,
                                                                 dataset="DMEM_6hr_screen_normal_pe_1",
                                                                 sigset="screen_large",
                                                                 method = "mygsea") {
  printCurrentFunction()
  npg <- 48
  if(do.read) {
    dir <- "../output/signature_conc_resp_summary/"
    file <- paste0(dir,"SIGNATURE_CR_screen_large_DMEM_6hr_screen_normal_pe_1_mygsea_0.05_conthits hits.RData")
    print(file)
    load(file)
    SCR <<- SIG_CONC_RESPONSE_ACTIVE
    file <- "../output/noiseHunter/httr_mcf7_ph1_metadata.xlsx"
    smap <- read.xlsx(file)
    SMAP <<- smap
  }
  scr <- SCR
  smap <- SMAP
  rownames(smap) <- smap$sample_id
  smap <- smap[is.element(smap$stype,"test sample"),]
  name.list <- c("dtxsid","casrn","name","sample_id","block_id","pg_id",
                 "nhit","nhit.10","nhit.1","nhit.0.1",
                 "n_reads","n_reads_mapd","mapd_frac","bad_probe_count","n_cov5","n_sig80","top10_prop","gini_coef",
                 "n_reads.sd","n_reads_mapd.sd","mapd_frac.sd","bad_probe_count.sd","n_cov5.sd","n_sig80.sd","top10_prop.sd","gini_coef.sd")
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

      res[k,"nhit"] <- nrow(temp)
      temp <- temp[temp$bmd<10,]
      res[k,"nhit.10"] <- nrow(temp)
      temp <- temp[temp$bmd<1,]
      res[k,"nhit.1"] <- nrow(temp)
      temp <- temp[temp$bmd<0.1,]
      res[k,"nhit.0.1"] <- nrow(temp)

      stemp <- smap.pg[is.element(smap.pg$chem_id,cid),]
      col.list <- names(res)[11:ncol(res)]
      col.list <- col.list[is.element(col.list,names(stemp))]

      for(col in col.list) {
        res[k,col] <- mean(stemp[,col])
        col.sd <- paste0(col,".sd")
        res[k,col.sd] <- sd(stemp[,col])
      }
    }
    result <- rbind(result,res)
    #browser()
  }
  file <- "../output/noiseHunter/noiseHunter.signatureHits.vs.rawCountsStats.xlsx"
  write.xlsx(result,file)
}
