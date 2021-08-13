#--------------------------------------------------------------------------------------
#' Generate chemicalwise boxplot of the BMD distributions by super_target
#'
#'   heparg2d_toxcast_pfas_pe1_normal
#'   mcf7_ph1_pe1_normal_block_123_allPG
#'   u2os_toxcast_pfas_pe1_normal
#'   PFAS_HepaRG
#'   PFAS_U2OS
#'   u2os_pilot_pe1_normal_null_pilot_lowconc
#'   u2os_toxcast_pfas_pe1_normal_refchems
#'   heparg2d_toxcast_pfas_pe1_normal_refchems
#'
#' MCF7_HTTr_Water_Samples
#'
#' After running this, run the following ...
#' superTargetPODplot
#' superTargetStats
#--------------------------------------------------------------------------------------
signatureAndGenePOD <- function(to.file=T,
                               do.load=T,
                               dataset="heparg2d_toxcast_pfas_pe1_normal_refchems",
                               sigcatalog="signatureDB_master_catalog 2021-04-24",
                               sigset="screen_large",
                               method="fc",
                               celltype="HepaRG",
                               hccut=0.95,
                               tccut=1.5,
                               minconc=0.0001,
                               maxconc=100) {
  printCurrentFunction(paste(dataset,sigset,method))

  if(to.file) {
    file = paste0("../output/pod/signatureAndGenePOD_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,".pdf")
    pdf(file=file,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,1),mar=c(4,15,6,2))

  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    SMAT <<- SIGNATURE_CR

    file <- paste0("../output/gene_conc_resp_summary/GENE_CR_",dataset,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    GMAT <<- GENE_CR
  }
  smat = SMAT
  gmat = GMAT

  smat = smat[smat$hitcall>hccut,]
  smat = smat[smat$top_over_cutoff>tccut,]
  gmat = gmat[gmat$hitcall>hccut,]
  gmat = gmat[gmat$top_over_cutoff>tccut,]

  name.list = c("sample_id","dtxsid","name","pod_gene","pod_signature","pod_super_target")
  sids = unique(c(smat$sample_id,gmat$sample_id))
  nsid = length(sids)
  res = as.data.frame(matrix(nrow=nsid,ncol=length(name.list)))
  names(res) = name.list
  res$sample_id = sids
  rownames(res) = sids
  res$pod_gene = -1
  res$pod_signature = -1
  res$pod_super_target = -1
  for(sid in sids) {
    temp = smat[is.element(smat$sample_id,sid),]
    res[sid,"dtxsid"] = temp[1,"dtxsid"]
    res[sid,"name"] = temp[1,"name"]
    if(nrow(temp)>2) {
      bmds = sort(temp$bmd)
      qb = quantile(bmds,probs=seq(0,1,0.05))
      res[sid,"pod_signature"] = qb[2]
    }

    st.list = unique(temp$super_target)
    nst = length(st.list)
    if(nst>0) {
      sres = as.data.frame(matrix(nrow=nst,ncol=2))
      sres[,2] = -1
      for(i in 1:nst) {
        st = st.list[i]
        sres[i,1] = st
        vals = temp[is.element(temp$super_target,st),"bmd"]
        if(length(vals)>2) sres[i,2] = median(vals)
      }
      vals = sres[,2]
      vals = vals[vals>0]
      if(length(vals)>0) res[sid,"pod_super_target"] = min(vals)

    }
    temp = gmat[is.element(gmat$sample_id,sid),]
    if(nrow(temp)>2) {
      bmds = sort(temp$bmd)
      qb = quantile(bmds,probs=seq(0,1,0.05))
      res[sid,"pod_gene"] = qb[2]
    }
  }

  x = res$name
  y = res$pod_gene
  x = x[y>0]
  y = y[y>0]
  main = paste(dataset,"\ngene")
  boxplot(y~x,main=main,cex.main=1.2,
          ylim=c(minconc,100),log="x",xlab="BMD (uM)",ylab="",
          horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0))
  for(i in -5:3) {
    x = 10**i
    lines(c(x,x),c(0,100))
  }
  x = res$name
  y = res$pod_signature
  x = x[y>0]
  y = y[y>0]
  main = paste(dataset,"\nsignature")
  boxplot(y~x,main=main,cex.main=1.2,
          ylim=c(minconc,100),log="x",xlab="BMD (uM)",ylab="",
          horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0))
  for(i in -5:3) {
    x = 10**i
    lines(c(x,x),c(0,100))
  }

  x = res$name
  y = res$pod_super_target
  x = x[y>0]
  y = y[y>0]
  main = paste(dataset,"\nsuper target")
  boxplot(y~x,main=main,cex.main=1.2,
          ylim=c(minconc,100),log="x",xlab="BMD (uM)",ylab="",
          horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0))
  for(i in -5:3) {
    x = 10**i
    lines(c(x,x),c(0,100))
  }
  if(!to.file) browser()
  else dev.off()

  file = paste0("../output/pod/signatureAndGenePOD_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,".xlsx")
  write.xlsx(res,file)


}

