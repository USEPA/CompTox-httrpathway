#--------------------------------------------------------------------------------------
#' Generate chemicalwise distribution plots by target group
#'
#'  heparg2d_toxcast_pfas_pe1_normal
#'  mcf7_ph1_pe1_normal_block_123
#'  u2os_toxcast_pfas_pe1_normal
#'  PFAS_HepaRG
#'  PFAS_U2OS
#'  u2os_pilot_pe1_normal_null_pilot_lowconc
#'  u2os_toxcast_pfas_pe1_normal_refchems
#'  heparg2d_toxcast_pfas_pe1_normal_refchems
#'
#' MCF7_HTTr_Water_Samples
#--------------------------------------------------------------------------------------
signatureDensityPlot <- function(to.file=F,
                                 do.load=F,
                                 dataset="u2os_toxcast_pfas_pe1_normal_refchems",
                                 sigset="screen_large",
                                 method="fc",
                                 celltype="U2OS",
                                 sigcatalog="signatureDB_master_catalog 2021-03-05",
                                 target_class="Immune",
                                 hccut=0.95,
                                 tccut=1.5) {
  printCurrentFunction(paste(dataset,sigset,method))

  file = paste0("../input/signatures/",sigcatalog,".xlsx")
  catalog = read.xlsx(file)
  catalog = catalog[catalog[,sigset]==1,]
  catalog = unique(catalog[,c("parent","super_target","target_class")])

  file = "../input/chemicals/PFAS synonyms.xlsx"
  synonyms = read.xlsx(file)
  synonyms = synonyms[!is.na(synonyms$nickname),]
  rownames(synonyms) = synonyms$dtxsid

  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    mat[is.na(mat$bmd),"bmd"] = 1000
    mat[is.na(mat$top_over_cutoff),"top_over_cutoff"] = 0
    mat = mat[mat$hitcall>hccut,]
    mat = mat[mat$top_over_cutoff>tccut,]
    MAT <<- mat
  }
  mat = MAT

  sig.list.all = catalog[,"parent"]
  nall = length(sig.list.all)
  sig.list.in = catalog[is.element(catalog$target_class,target_class),"parent"]
  sig.list.out = catalog[!is.element(catalog$target_class,target_class),"parent"]

  temp = mat[is.element(mat$name,"Cyclosporin A"),]
  temp = temp[is.element(temp$signature,sig.list.in),]
  if(length(temp$signature)>0) sig.list.in = unique(temp$signature)

  if(to.file) {
    fname <- paste0("../output/signature_density_plots/signature_density_plot",target_class,"_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,5,5,2))

  res= unique(mat[,c("name","dtxsid","sample_id")])
  rownames(res) = res$sample_id
  res$burstlimit=1000
  chem.list = sort(unique(mat$name))
  for(chem in chem.list) {
    cat(chem,"\n")
    tempchem = mat[is.element(mat$name,chem),]
    dtxsid = tempchem[1,"dtxsid"]
    if(is.element(dtxsid,synonyms$dtxsid)) {
      nn = synonyms[dtxsid,"nickname"]
      chem = paste(nn,":",chem)
    }

    bw = 0.05
    sid.list = unique(tempchem$sample_id)
    for(sid in sid.list) {
      temp = tempchem[is.element(tempchem$sample_id,sid),]
      bmd.all =  log10(temp[,"bmd"])
      bmd.in =  log10(temp[is.element(temp$signature,sig.list.in),"bmd"])
      bmd.out = log10(temp[is.element(temp$signature,sig.list.out),"bmd"])
      if(length(bmd.in)>2 && length(bmd.out)>2) {
        #########################################################################
        x = density(bmd.all,bw=bw)
        peak = which.max(x$y)
        ymax = 2*max(x$y)
        peak.loc = x$x[peak]
        plot(x,xlim=c(-2,3),ylim=c(0,ymax),xlab="log(bmd uM)",ylab="",cex.lab=1.2,cex.axis=1.2,
             main=paste(chem,"\n",dtxsid,sid),type="n")
        left = peak.loc-1
        right = 3
        rect(left,0,right,ymax,col="lightgray")
        lines(density(bmd.all,bw=bw),col="black")
        lines(c(peak.loc,peak.loc),c(0,10))
        frac = 100*length(bmd.all)/nall
        text(-2,0.9*ymax,paste("Hit:",format(frac,digits=2),"%"),pos=4,cex=1.5)

        #########################################################################
        x = density(bmd.out,bw=bw)
        peak = which.max(x$y)
        peak.loc = x$x[peak]
        plot(x,xlim=c(-2,3),ylim=c(0,ymax),xlab="log(bmd uM)",ylab="",cex.lab=1.2,cex.axis=1.2,
             main=paste(chem,"\n",dtxsid,sid),type="n")
        left = peak.loc-1
        right = 3
        rect(left,0,right,ymax,col="lightgray")
        lines(density(bmd.out,bw=bw),col="black")
        lines(density(bmd.in,bw=bw),col="red")
        lines(c(peak.loc,peak.loc),c(0,10))
        res[sid,"burstlimit"] = 10**left
        text(-2,0.9*ymax,paste("Hit:",format(frac,digits=2),"%"),pos=4,cex=1.5)
        #########################################################################
        if(!to.file) browser()
      }
    }
  }
  file = paste0("../output/signature_density_plots/signature_density_plot",target_class,"_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,"_burstlimit.xlsx")
  write.xlsx(res,file)

  if(to.file) dev.off()
}

