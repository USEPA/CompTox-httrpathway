#--------------------------------------------------------------------------------------
#' Generate distributions of the potency for on and off target activity
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#' u2os_pilot_pe1_normal_null_pilot_lowconc
#' u2os_toxcast_pfas_pe1_normal_refchems
#' heparg2d_toxcast_pfas_pe1_normal_refchems
#'
#--------------------------------------------------------------------------------------
superTargetOnOffTarget <- function(to.file=F,
                               do.load=F,
                               dataset="mcf7_ph1_pe1_normal_block_123",
                               sigset="screen_large",
                               method="fc",
                               celltype="MCF7",
                               hccut=0.95,
                               tccut=1.5,
                               minconc=0.0001,
                               maxconc=100) {
  printCurrentFunction(paste(dataset,sigset,method))

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
  if(to.file) {
    fname <- paste0("../output/super_target_boxplot/",celltype,"/super_target_on_off_target_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(5,4,4,2))
  mat = MAT
  file = "../input/chemicals/refchem_super_target_map.xlsx"
  stdb = read.xlsx(file)
  file = paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,"_stats.xlsx")
  stats = read.xlsx(file)
  stats = stats[stats$TP>=4,]
  stats = stats[stats$Sens>=0.5,]
  st.list = stats$rowname
  st.list = sort(st.list[is.element(st.list,mat$super_target)])

  #st.list =st.list[1:6]
  x.on = NULL
  y.on = NULL
  x.off = NULL
  y.off = NULL
  for(st in st.list) {
    cat(st,"\n")
    temp = mat[is.element(mat$super_target,st),]
    dtxsid.list = unique(temp$dtxsid)
    for(dtxsid in dtxsid.list) {
      dtemp = temp[is.element(temp$dtxsid,dtxsid),]
      y = dtemp$bmd
      x = y
      x[] = st
      stemp = stdb[is.element(stdb$dtxsid,dtxsid),]
      stl = str_split(stemp[1,"stlist"],"\\|")[[1]]
      if(is.element(st,stl)) {
        x.on = c(x.on,x)
        y.on = c(y.on,y)
      }
      else {
        x.off = c(x.off,x)
        y.off = c(y.off,y)
      }
    }
  }

  for(st in st.list) {
    val.on = log10(y.on[is.element(x.on,st)])
    val.off = log10(y.off[is.element(x.off,st)])
    p = wilcox.test(val.on,val.off)$p.value
    plot(density(val.off),xlim=c(-3,2),xlab="log(bmd (uM))",cex.lab=1.3,cex.axis=1.3,main=paste(st,"\np=",format(p,digits=2)),ylim=c(0,1.5))
    lines(density(val.on),col="red")
    if(!to.file) browser()
  }

  if(to.file) dev.off()
}

