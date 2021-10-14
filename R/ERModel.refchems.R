#--------------------------------------------------------------------------------------
#' Analyze the ER model data for reference chemicals
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

#' After running this, run the following ...
#' superTargetPODplot
#' superTargetStats
# dataset="mcf7_ph1_pe1_normal_block_123_allPG",
# sigcatalog="signatureDB_master_catalog 2021-08-27",
# sigset="screen_large",
#
# sigcatalog="signatureDB_master_catalog ER",sigset="estrogen",
#
#-------------------------------------------------------------------------------
ERModel.refchems <- function(to.file=T,
                             dataset="mcf7_ph1_pe1_normal_block_123_allPG",
                             sigcatalog="signatureDB_master_catalog 2021-08-27",
                             sigset="screen_large",
                             method="gsea",
                             celltype="MCF7",
                             hccut=0.9,
                             tccut=1,
                             minhit=10) {
  printCurrentFunction(paste(dataset,sigset,method))
  file = paste0("../ERModel/ER_chems all ",dataset," ",sigset," ",hccut," ",minhit,".xlsx")
  print(file)
  erchems = read.xlsx(file)

  #erchems[is.element(erchems$hts.mode,"inactive"),"hts.pod.agonist"] = 3
  #erchems[is.element(erchems$hts.mode,"inactive"),"hts.pod.antagonist"] = 3

  #erchems[erchems$nhit<2,"httr.mode"] = "inactive"
  #erchems[erchems$nhit<2,"httr.pod"] = 3
  if(to.file) {
    fname <- paste0("../ERModel/ERmodel refchems ",dataset," ",sigset," ",hccut," ",minhit,".pdf")
    pdf(file=fname,width=8,height=8,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }

  file = paste0("../ERModel/ER_chems burst ",dataset," ",sigset," ",hccut," ",minhit,".xlsx")
  burst = read.xlsx(file)
  rownames(burst) = burst$dtxsid

  ################################################################################################
  # Plot ER model AUC vs HTTr pod
  ################################################################################################
  par(mfrow=c(1,1),mar=c(4,4,4,2))
  temp = erchems[!is.na(erchems$auc.agonist),]
  temp = temp[!is.na(temp$auc.antagonist),]
  x = temp$auc.agonist
  for(i in 1:nrow(temp)) if(temp[i,"auc.antagonist"]>temp[i,"auc.agonist"]) x[i] = temp[i,"auc.antagonist"]
  y = temp$httr.pod
  mod = lm(y~x)
  smod = summary(mod)
  intercept = smod$coefficients[1,1]
  slope = smod$coefficients[2,1]
  p = smod$coefficients[2,4]
  plot(y~x,xlab="max(AUC(agonist),AUC(antagonist))",ylab="Mean(log(bmd uM))",cex.lab=1.2,cex.axis=1.2,ylim=c(-3,3),main=paste("p=",format(p,digits=2)))
  for(i in 1:nrow(temp)) {
    ag = temp[i,"auc.agonist"]
    antag = temp[i,"auc.antagonist"]
    col = "blue"
    if(antag>ag) col = "red"
    points(x[i],y[i],pch=21,bg=col)
  }
  lines(c(0,1.5),c(intercept,intercept+slope*1.5))
  #if(!to.file) browser()

  ################################################################################################
  # Plot ER pod AUC vs HTTr pod
  ################################################################################################
  par(mfrow=c(1,1),mar=c(4,4,4,2))
  x = erchems$hts.pod.agonist
  y = erchems$httr.pod

  xymin = -5
  plot(y~x,xlab="mean(log(HTS POD uM))",ylab="Mean(log(HTTr BMD uM))",cex.lab=1.2,cex.axis=1.2,xlim=c(xymin,3),ylim=c(xymin,3),
       main=paste0(dataset," ",sigset," ",hccut," ",minhit),type="n")

  lista = NULL
  listb = NULL
  for(i in 1:nrow(erchems)) {
    x = min(erchems[i,"hts.pod.agonist.10"],erchems[i,"hts.pod.antagonist.10"])
    y = erchems[i,"httr.pod"]
    if(!is.na(x) && !is.na(y)) {
      dtxsid = erchems[i,"dtxsid"]
      name = erchems[i,"name"]
      cat(dtxsid,name,x,y,"\n")
      if(dtxsid=="DTXSID8022325") name = "HPTE"
      col = "white"
      pch = 21
      mode1 = erchems[i,"httr.mode"]
      mode2 = erchems[i,"hts.mode"]
      if(mode1=="agonist") {
        if(is.na(mode2)) col = "blue"
        else {
          if(mode2=="agonist") col = "blue"
          else if(mode2=="antagonist") col = "cyan"
          else col = "gray"
        }
      }

      else if(mode1=="antagonist") {
        if(is.na(mode2)) col = "red"
        else {
          if(mode2=="antagonist") col = "red"
          else if(mode2=="agonist") col = "orange"
          else col = "gray"
        }
      }

      else if(mode1=="inactive") {
        if(is.na(mode2)) col = "white"
        else {
          if(mode2=="inactive") col = "white"
          else col = "green"
        }
      }

      if(is.element(dtxsid,burst$dtxsid)) {
        interference = burst[dtxsid,"interference"]
        if(interference==1) pch = 25
      }
      if(x==3) {
        if(col=="white") x = 2.7+rnorm(1,0,0.05)
        else x = x+rnorm(1,0,0.05)
      }
      if(y==3) {
        if(col=="white") y = 2.7+rnorm(1,0,0.05)
        else y = y+rnorm(1,0,0.05)
      }

      points(x,y,pch=pch,bg=col)
      if(col!="white") {
        #if((y>x+1) && x<1 && y<2 && col!="gray") text(x,y,name,pos=2,cex=0.7)
        #if((y<x-1) && x<1 && y<2 && col!="gray") text(x,y,name,pos=4,cex=0.7)
        #if(x>1 && y< -0.5) text(x,y,name,pos=2,cex=0.7)
      }
      if(x<2 && y>2 && col!="white") lista = c(lista,name)
      if(x>2 && y<0 && col!="white") listb = c(listb,name)
    }
  }
  x = xymin
  y = 1.9
  dy=0.2
  for(i in 1:length(lista)) {
    text(x,y,lista[i],pos=4,cex=0.7)
    y = y-dy
  }

  x = xymin
  y = 0
  dy=0.2
  for(i in 1:length(listb)) {
    text(x,y,listb[i],pos=4,cex=0.7)
    y = y-dy
  }

  lines(c(xymin,3),c(xymin,3))
  lines(c(xymin,3),c(xymin+1,4),col="gray")
  lines(c(xymin,3),c(xymin-1,2),col="gray")

  lines(c(-10,10),c(2,2),col="black")
  lines(c(2,2),c(-10,10),col="black")

  do.legend=T
  if(do.legend) {
    x = -1.5
    y = -3.5
    dy = 0.2
    points(x,y,pch=21,bg="blue")
    text(x,y,"Both Agonist",pos=4)
    y = y-dy
    points(x,y,pch=21,bg="red")
    text(x,y,"Both Antagonist",pos=4)
    y = y-dy
    points(x,y,pch=21,bg="white")
    text(x,y,"Both Inactive",pos=4)
    y = y-dy
    points(x,y,pch=21,bg="cyan")
    text(x,y,"HTTr Agonist / HTS Antagonist",pos=4)
    y = y-dy
    points(x,y,pch=21,bg="green")
    text(x,y,"HTTr Inactive, HTS Active",pos=4)
    y = y-dy
    points(x,y,pch=21,bg="gray")
    text(x,y,"HTTr Active, HTS Inactive",pos=4)
    y = y-dy
    points(x,y,pch=25,bg="gray")
    text(x,y,"HTTr Assay Interference",pos=4)
  }

  text(-4.6,-5,"A",cex=2)
  text(-4.6,2.3,"B",cex=2)
  text(3,-5,"C",cex=2)
  text(3,2.3,"D",cex=2)

  if(!to.file) browser()

  ################################################################################################
  # Plot Referecne chemical stats
  ################################################################################################
  par(mfrow=c(3,2),mar=c(4,10,4,2))

  col.list = c("refchem.invitro.agonist","refchem.invitro.antagonist","refchem.invivo.agonist")
  header.list = c("In Vitro (agonist)","In Vitro (antagonist)","In Vivo (agonist)")
  for(j in 1:3) {
    col = col.list[j]
    header = header.list[j]
    temp = erchems[!is.na(erchems[,col]),]
    if(col=="refchem.invitro.agonist") temp = temp[is.element(temp$httr.mode,"agonist"),]
    if(col=="refchem.invitro.antagonist") temp = temp[is.element(temp$httr.mode,"antagonist"),]
    if(col=="refchem.invivo.agonist") temp = temp[is.element(temp$httr.mode,"agonist"),]
    x = NULL
    y = NULL
    for(i in 1:nrow(temp)) {
      #if(col=="refchem.invitro.agonist") label = paste(temp[i,col],":",temp[i,"mode"])
      #else
      label = temp[i,col]
      value = 10**temp[i,"httr.pod"]
      x = c(x,label)
      y = c(y,value)

    }
    if(col=="refchem.invitro.agonist") {
      x[is.element(x,"Inactive")]  = "    Inactive"
      x[is.element(x,"Very Weak")] = "   Very Weak"
      x[is.element(x,"Weak")]      = "  Weak"
      x[is.element(x,"Moderate")]  = " Moderate"
      x[is.element(x,"Strong")]    = "Strong"
    }
    boxplot(y~x,main=header,cex.main=0.8,
            ylim=c(1e-3,1e3),log="x",xlab="BMD (uM)",ylab="",
            horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0))
    for(x in c(1e-4,1e-3,1e-2,1e-1,1,10,100,1000)) lines(c(x,x),c(0,100),col="gray")
  }
  if(!to.file) browser()

  if(!to.file) browser()
  else dev.off()
}
