#--------------------------------------------------------------------------------------
#' Generate chemicalwise boxplot of the BMD distributions by super_target
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#'
#--------------------------------------------------------------------------------------
plateGroupBoxplot <- function(to.file=F,
                               do.load=F,
                               dataset="mcf7_ph1_pe1_normal_block_123",
                               sigset="screen_large",
                               method="fc",
                               celltype="MCF7",
                               hccut=0.9,
                               tccut=1.5) {
  printCurrentFunction(paste(dataset,sigset,method))
  file = paste0("../input/fcdata/spid.to.pg.map ",dataset,".xlsx")
  map = read.xlsx(file)

  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    MAT <<- mat
  }
  mat = MAT
  mat[is.na(mat$bmd),"bmd"] = 1000
  mat[is.na(mat$top_over_cutoff),"top_over_cutoff"] = 0
  mat[mat$hitcall<hccut,"bmd"] = 1000
  mat[mat$hitcall<hccut,"hitcall"] = 0
  mat[mat$top_over_cutoff<tccut,"bmd"] = 1000
  mat[mat$top_over_cutoff<tccut,"hitcall"] = 0
  mat = mat[mat$hitcall>0,]
  mat = mat[mat$bmd<1,]
  print(nrow(MAT))
  print(nrow(mat))
  x = NULL
  y = NULL
  pg.list = unique(map$pg_id)
  for(pg_id in pg.list) {
    spid.list = map[is.element(map$pg_id,pg_id),"spid"]
    yy = mat[is.element(mat$sample_id,spid.list),"bmd"]
    xx = yy
    xx[] = pg_id
    x = c(x,xx)
    y = c(y,yy)
  }

  if(to.file) {
    fname <- paste0("../output/super_target_boxplot/",celltype,"/plateGroupBoxplot",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(1,1),mar=c(4,5,4,2))
  boxplot(y~x,main="Low Conc BMD by Plate Group",
          ylim=c(0.00001,1),log="x",xlab="BMD (uM)",ylab="",
          horizontal=T,las=1,par(cex.lab=0.8,cex.axis=0.8))
  if(!to.file) browser()
  else dev.off()

}

