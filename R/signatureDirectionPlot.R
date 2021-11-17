#--------------------------------------------------------------------------------------
#' Plot the cummulative distribution functions of the up and down direction signatures
#'
#' @param to.file If TRUE, send the plots to a file
#' @param do.load If TRUE, load hte large HTTr data set into memory
#' @param dataset Name of the HTTr data set
#' @param sigcatalog Name of the signature catalog to use
#' @param sigset Name of the signature set
#' @param method Scoring method
#' @param celltype Name of the cell type
#' @param hccut Exclude rows in the data set with hitcall less than this value
#' @param tccut Exclude rows in the data set with top_over_cutoff less than this value
#' @param cutoff The minimum number of signatures hat have to be active in a super
#' target for the super target to be considered active. Default is 5
#' @param minconc Minimum concentration used in the plots
#' @param maxconc Maximum concentration used in the plots
#'
#' After running this, run the following ...
#' superTargetPODplot
#' superTargetStats
#--------------------------------------------------------------------------------------
signatureDirectionPlot <- function(to.file=T,
                                   do.load=F,
                                   dataset="MCF7_pilot_PRF_6hr_pilot_normal_pe_1",
                                   sigset="screen_large",
                                   method="gsea",
                                   celltype="MCF7",
                                   hccut=0.9,
                                   tccut=1) {
  printCurrentFunction(paste(dataset,sigset,method))

  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    mat = mat[!is.na(mat$bmd),]
    mat = mat[!is.na(mat$top_over_cutoff),]
    mat[mat$hitcall<hccut,"bmd"] = 1000
    mat[mat$hitcall<hccut,"hitcall"] = 0
    mat[mat$top_over_cutoff<tccut,"bmd"] = 1000
    mat[mat$top_over_cutoff<tccut,"hitcall"] = 0
    mat = mat[mat$hitcall>0,]
    MAT <<- mat
  }
  mat = MAT

  if(to.file) {
    fname = paste0("../output/signature_direction_plot/signature_direction_plot ",celltype," ",dataset," ",sigset," ",method," ",hccut," ",tccut,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(4,3),mar=c(3,3,3,2))
  dlist = unique(mat$dtxsid)
  xgrid = seq(from=-3,to=3,by=0.1)
  xgrid = 10**xgrid

  up = xgrid
  up[] = 0
  dn = up
  for(i in 1:length(xgrid)) {
    conc = xgrid[i]
    temp1 = mat[mat$bmd<conc,"top"]
    up[i] = length(temp1[temp1>0])
    dn[i] = length(temp1[temp1<0])
  }
  mval = max(c(up,dn))
  up = up/mval
  dn = dn/mval
  plot(up~xgrid,log="x",xlim=c(1e-2,1e2),ylim=c(0,1),main="All Chemicals",cex.axis=1.2,cex.lab=1.2,
       type="l",col="red",xlab="Conc",ylab="",cex.main=0.8)
  lines(dn~xgrid,col="black")



  for(dtxsid in dlist) {
    temp = mat[is.element(mat$dtxsid,dtxsid),]
    name = temp[1,"name"]
    up = xgrid
    up[] = 0
    dn = up
    for(i in 1:length(xgrid)) {
      conc = xgrid[i]
      temp1 = temp[temp$bmd<conc,"top"]
      up[i] = length(temp1[temp1>0])
      dn[i] = length(temp1[temp1<0])
    }
    mval = max(c(up,dn))
    up = up/mval
    dn = dn/mval
    plot(up~xgrid,log="x",xlim=c(1e-2,1e2),ylim=c(0,1),main=name,cex.axis=1.2,cex.lab=1.2,
         type="l",col="red",xlab="Conc",ylab="",cex.main=0.8)
    lines(dn~xgrid,col="black")
    if(!to.file) browser()
  }
  if(to.file) dev.off()

}

