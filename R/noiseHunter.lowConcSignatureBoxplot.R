#--------------------------------------------------------------------------------------
#' Look for trends in the raw data by plate within plate group
#'
#' @param dir The directory from which to read all of the raw filesatalog file
#--------------------------------------------------------------------------------------
noiseHunter.lowConcSignatureBoxplot <- function(to.file=F,
                                                do.read=F,
                                                dataset="DMEM_6hr_screen_normal_pe_1",
                                                sigset="screen_large",
                                                method = "mygsea") {
  printCurrentFunction()
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
  if(to.file) {
    fname <- paste0("../output/noiseHunter/noiseHunter.lowConcSignatureBoxplot.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(1,1),mar=c(5,4,4,3))

  scr <- SCR
  smap <- SMAP

  pg.list <- NULL
  bmd.list <- NULL
  npg <- 48
  #npg <- 4

  for(pg in 1:npg) {
    cat(pg,"\n")
    cid.list <- smap[smap$pg_id==pg,"chem_id"]
    temp <- scr[is.element(scr$sample_id,cid.list),]
    #temp <- temp[1:1000,]
    bmd <- log10(temp$bmd)
    pgs <- bmd
    pgs[] <- pg
    pg.list <- c(pg.list,pgs)
    bmd.list <- c(bmd.list,bmd)
  }
  boxplot(bmd.list~pg.list,horizontal=T,cex.lab=1.2,cex.axis=1,
          ylab="Plate group",xlab="log10(BMD)",las=2,outcex=0.1,outpch=21,bg="gray",fg="gray",col="gray")

  if(to.file) dev.off()
  else browser()
}
