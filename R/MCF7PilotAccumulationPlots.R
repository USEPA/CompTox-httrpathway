#' Build accumulation plots for ER and HMGCR
#'
MCF7PilotAccumulationPlots <- function(to.file=F,
                              method="gsea",
                              celltype="MCF7",
                              sigset="screen_large",
                              hccut=0.9,
                              tccut=1) {
  printCurrentFunction()
  dir = "../output/mcf7_pilot/"
  if(to.file) {
    fname = paste0(dir,"MCF7PilotAccumulationPlots",method,"_",hccut,"_",tccut,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,5,6,2))

  ######################################################################################################
  # read the httrpathway data
  ######################################################################################################
  ####################################################################################
  # compare the PODs to the ER model values
  ####################################################################################
  file = paste0(dir,"ER_chems all mcf7_ph1_pe1_normal_block_123_allPG screen_large 0.9 10.xlsx")
  ermodel = read.xlsx(file)
  erchems = c(
    "Fulvestrant",
    "4-Hydroxytamoxifen",
    "Clomiphene citrate (1:1)",
    "Bisphenol A",
    "Bisphenol B",
    "4-Nonylphenol, branched",
    "4-Cumylphenol"
  )
  ermodel = ermodel[is.element(ermodel$name,erchems),]
  rownames(ermodel) = ermodel$name
  ermodel = ermodel[erchems,]
  rownames(ermodel) = ermodel$dtxsid
  nchem = nrow(ermodel)
  dlist = ermodel$dtxsid


  nset = 6
  dataset.list = c(
    "MCF7_pilot_DMEM_6hr_pilot_normal_pe_1",
    "MCF7_pilot_DMEM_12hr_pilot_normal_pe_1",
    "MCF7_pilot_DMEM_24hr_pilot_normal_pe_1",
    "MCF7_pilot_PRF_6hr_pilot_normal_pe_1",
    "MCF7_pilot_PRF_12hr_pilot_normal_pe_1",
    "MCF7_pilot_PRF_24hr_pilot_normal_pe_1"
  )
  media.list = c("DMEM","DMEM","DMEM","PRF","PRF","PRF")
  time.list = c(6,12,24,6,12,24)
  res = NULL
  for(i in 1:nset) {
    par(mfrow=c(4,3),mar=c(4,5,4,2))
    dataset = dataset.list[i]
    media = media.list[i]
    time = time.list[i]
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    res = SIGNATURE_CR
    res = res[res$hitcall>hccut,]
    res = res[res$top_over_cutoff>tccut,]
    for(dtxsid in dlist) {
      temp1 = res[is.element(res$dtxsid,dtxsid),]
      name = temp1[1,"name"]
      cat(name,"\n")
      st = "Estrogen"
      pod.in = temp1[is.element(temp1$super_target,st),"bmd"]
      pod.out = temp1[!is.element(temp1$super_target,st),"bmd"]
      if(length(pod.in)>0) {
        x = seq(from=-3,to=2,by=0.2)
        x = 10**x
        y1 = x
        y1[] = 0
        y2 = y1
        for(j in 1:length(x)) {
          conc = x[j]
          y1[j] = length(pod.in[pod.in<conc])/length(pod.in)
          y2[j] = length(pod.out[pod.out<conc])/length(pod.out)
        }
        plot(y2~x,main=paste(name,"\n",media," ",time,"hr"),type="l",log="x",xlim=c(1e-3,100),ylim=c(0,1),xlab="Conc (uM)",ylab="Cumulative Dist",cex.lab=1.2,cex.axis=1.2,col="black")
        lines(y1~x,col="red")

        pod.hts = min(ermodel[dtxsid,"hts.pod.agonist"],ermodel[dtxsid,"hts.pod.antagonist"])
        pod.hts = 10**pod.hts
        lines(c(pod.hts,pod.hts),c(0,1),col="red",lwd=2)
       }
    }
    for(dtxsid in c("DTXSID0023581","DTXSID5020784")) {
      temp1 = res[is.element(res$dtxsid,dtxsid),]
      name = temp1[1,"name"]
      cat(name,"\n")
      temp1 = temp1[temp1$hitcall>hccut,]
      temp1 = temp1[temp1$top_over_cutoff>tccut,]
      st = c("HMGCR","Cholesterol")
      pod.in = temp1[is.element(temp1$super_target,st),"bmd"]
      pod.out = temp1[!is.element(temp1$super_target,st),"bmd"]
      if(length(pod.in)>0) {
        x = seq(from=-3,to=2,by=0.2)
        x = 10**x
        y1 = x
        y1[] = 0
        y2 = y1
        for(j in 1:length(x)) {
          conc = x[j]
          y1[j] = length(pod.in[pod.in<conc])/length(pod.in)
          y2[j] = length(pod.out[pod.out<conc])/length(pod.out)
        }
        plot(y2~x,main=paste(name,"\n",media," ",time,"hr"),type="l",log="x",xlim=c(1e-3,100),ylim=c(0,1),xlab="Conc (uM)",ylab="Cumulative Dist",cex.lab=1.2,cex.axis=1.2,col="black")
        lines(y1~x,col="red")
        if(name=="Simvastatin") pod.hts = 11.2
        if(name=="Lovastatin") pod.hts = 8.4
        lines(c(pod.hts,pod.hts),c(0,1),col="red",lwd=2)
      }
    }
    if(!to.file) browser()
  }
  if(to.file) dev.off()
}
