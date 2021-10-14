#--------------------------------------------------------------------------------------
#' Look at concentration-slides of gene CR data to understand where burst starts
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
# dataset="mcf7_ph1_pe1_normal_block_123_allPG",
# sigcatalog="signatureDB_master_catalog 2021-08-27",
# sigset="screen_large",
#
#
# sigcatalog="signatureDB_master_catalog ER",sigset="estrogen"
#
#
#-------------------------------------------------------------------------------
signatureSlice <- function(to.file=F,
                           do.load=F,
                           dataset="mcf7_ph1_pe1_normal_block_123_allPG",
                           sigcatalog="signatureDB_master_catalog 2021-09-29",
                           sigset="screen_large",
                           method="gsea",
                           celltype="MCF7",
                           hccut=0.9,
                           minhit=10,
                           tccut = 0.9) {
  printCurrentFunction(paste(dataset))
  if(to.file) {
    fname <- paste0("../ERModel/signatureSlice ",dataset," ",sigset," ",hccut," ",minhit,".pdf")
    pdf(file=fname,width=8,height=8,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,5,5,2))

  ############################################################################################################
  # read the data
  ############################################################################################################
  file = paste0("../ERModel/ER_chems all ",dataset," ",sigset," ",hccut," ",minhit,".xlsx")
  erchems = read.xlsx(file)
  erchems = erchems[erchems$candidate==1,]
  rownames(erchems) = erchems$dtxsid
  if(do.load) {
    cat("Load the data\n")
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    temp = SIGNATURE_CR
    dlist = erchems$dtxsid
    temp = temp[is.element(temp$dtxsid,dlist),]
    SIG_CR <<- temp
  }
  print(dim(SIG_CR))

  file = "../output/cytotox/2021-08-07_MCF7_Screen_Viability_Data.RData"
  load(file=file)
  cytotox = viability.data

  burst = erchems[,c("dtxsid","name","httr.pod")]
  burst$peak.up = NA
  burst$shoulder.up = NA
  burst$diff.up = NA
  burst$peak.dn = NA
  burst$shoulder.dn = NA
  burst$diff.dn = NA
  burst$interference = 0

  sigmat = SIG_CR
  nsig = length(unique(SIG_CR$signature))
  sigmat$scaled_top = sigmat$top_over_cutoff * sign(sigmat$top)
  sigmat = sigmat[sigmat$top_over_cutoff>tccut,]
  sigmat = sigmat[!is.na(sigmat$bmd),]
  cnames = erchems$name
  breaks = seq(from=-4,to=3,by=0.25)
  sgenes.up = NULL
  sgenes.dn = NULL
  sgenes.low.up = NULL
  sgenes.low.dn = NULL
  for(cname in cnames) {
    temp = sigmat[is.element(sigmat$name,cname),]
    temp$logbmd = log10(temp$bmd)
    dtxsid = temp[1,"dtxsid"]
    pod = erchems[dtxsid,"httr.pod"]
    temp.up = temp[temp$top>0,]
    vals = temp.up$logbmd
    vals = vals[vals>breaks[1]]
    vals = vals[vals<breaks[length(breaks)]]
    h = hist(vals,breaks=breaks,plot=F)
    y = h$density
    x = breaks[1:length(y)]
    y = y/max(y)
    itop1 = which(y==1)[1]
    itop2 = itop1 - 2
    itop3 = itop2 - 2
    temp1.up = temp.up[temp.up$logbmd<=x[itop1],]
    temp1.up = temp1.up[temp1.up$logbmd>=x[itop2],]
    burst[dtxsid,"peak.up"] = x[itop1]
    burst[dtxsid,"shoulder.up"] = x[which.max(y)[1]]
    burst[dtxsid,"diff.up"] = burst[dtxsid,"peak.up"] - burst[dtxsid,"httr.pod"]


    temp2.up = temp.up[temp.up$logbmd<x[itop3],]
    sgenes.low.up = c(sgenes.low.up,unique(temp2.up$gene))

    if(nrow(temp1.up)>=5) {
      if(pod<x[itop2]-1) {
        sgenes.up = c(sgenes.up,unique(temp1.up$gene))
      }
      cat(cname,"up",nrow(temp1.up),"\n")

      plot(y~x,main=paste0(cname,"\n",nrow(temp1.up)),cex.lab=1.2,cex.axis=1.2,xlab="BMD (uM)",type="l",xlim=c(-3,3),ylim=c(0,1.2))
      lines(c(-5,5),c(1,1))
      lines(c(pod,pod),c(0,100),lwd=2,col="red")
      lines(c(x[itop1],x[itop1]),c(0,100))
      lines(c(x[itop2],x[itop2]),c(0,100))
      lines(c(x[itop3],x[itop3]),c(0,100),col="gray")

      temp.dn = temp[temp$top<0,]
      if(nrow(temp.dn)>0) {
        vals = temp.dn$logbmd
        vals = vals[vals>breaks[1]]
        vals = vals[vals<breaks[length(breaks)]]
        h = hist(vals,breaks=breaks,plot=F)
        y = h$density
        x = breaks[1:length(y)]
        y = y/max(y)
        itop1 = which(y==1)[1]
        itop2 = itop1 - 2

        if(is.na(itop1)) browser()
        burst[dtxsid,"peak.dn"] = x[itop1]
        burst[dtxsid,"shoulder.dn"] = x[which.max(y)[1]]
        burst[dtxsid,"diff.dn"] = burst[dtxsid,"peak.up"] - burst[dtxsid,"httr.pod"]

        temp1.dn = temp.dn[temp.dn$logbmd<=x[itop1],]
        temp1.dn = temp1.dn[temp1.dn$logbmd>=x[itop2],]
        if(pod<x[itop2]-1) sgenes.dn = c(sgenes.dn,unique(temp1.dn$gene))

        temp2.dn = temp.dn[temp.dn$logbmd<x[itop3],]
        sgenes.low.dn = c(sgenes.low.dn,unique(temp2.dn$gene))

        cat(cname,"dn",nrow(temp1.dn),"\n")
        lines(y~x,col="cyan")

        z = cytotox[is.element(cytotox$dtxsid,dtxsid),]
        x1 = log10(z$conc)
        y1 = z$norm_cell_count / 100
        points(y1~x1,pch=21,bg="red")

        z = cytotox[is.element(cytotox$dtxsid,dtxsid),]
        x1 = log10(z$conc)
        y1 = z$percent_resp_pi / 100
        points(y1~x1,pch=21,bg="cyan")

        z = cytotox[is.element(cytotox$dtxsid,dtxsid),]
        x1 = log10(z$conc)
        y1 = z$percent_resp_casp / 100
        points(y1~x1,pch=21,bg="gray")


        if(burst[dtxsid,"diff.up"]<1 || burst[dtxsid,"diff.dn"]<1) burst[dtxsid,"interference"] = 1
        if(burst[dtxsid,"interference"]==1) text(-3,0,"X",pos=4,cex=2)
      }
      if(!to.file) browser()
    }
  }
  genetable = table(sgenes.up)
  genetable = as.data.frame(sort(genetable,decreasing=T))
  file = paste0("../ERModel/ER_chems sigslice.up ",dataset," ",sigset," ",hccut," ",minhit,".xlsx")
  write.xlsx(genetable,file)
  genetable = table(sgenes.dn)
  genetable = as.data.frame(sort(genetable,decreasing=T))
  file = paste0("../ERModel/ER_chems sigslice.dn ",dataset," ",sigset," ",hccut," ",minhit,".xlsx")
  write.xlsx(genetable,file)

  genetable = table(sgenes.low.up)
  genetable = as.data.frame(sort(genetable,decreasing=T))
  file = paste0("../ERModel/ER_chems sigslice.low.up ",dataset," ",sigset," ",hccut," ",minhit,".xlsx")
  write.xlsx(genetable,file)

  genetable = table(sgenes.low.dn)
  genetable = as.data.frame(sort(genetable,decreasing=T))
  file = paste0("../ERModel/ER_chems sigslice.low.dn ",dataset," ",sigset," ",hccut," ",minhit,".xlsx")
  write.xlsx(genetable,file)

  file = paste0("../ERModel/ER_chems burst ",dataset," ",sigset," ",hccut," ",minhit,".xlsx")
  write.xlsx(burst,file)

  if(!to.file) browser()
  else dev.off()
}
