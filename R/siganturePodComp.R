#--------------------------------------------------------------------------------------
#'
#' Compare PODs from different conditions
#'
#' @param bmd.mode percent or abs
#' @param to.file If TRUE, write plots to a file
#--------------------------------------------------------------------------------------
signaturePodComp <- function(to.file=F,
                             sigset="screen_large",
                             method="fc") {
  printCurrentFunction()
  #' heparg2d_toxcast_pfas_pe1_normal
  #' mcf7_ph1_pe1_normal_good_pg
  #' u2os_toxcast_pfas_pe1_normal
  file = "../input/chemicals/httr_chemical_annotations 2020-09-14.xlsx"
  cann = read.xlsx(file)
  cann[is.na(cann$use_class),"use_class"] = "Other"
  if(to.file) {
    fname <- paste0("../output/super_target_boxplot/signaturePodComp_",sigset,"_",method,".pdf")
    pdf(file=fname,width=7,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,4,4))
  par(mfrow=c(2,1),mar=c(4,4,4,4))

  col.list = c("red","orange","yellow","green","blue","violet","gray","white","black","brown")
  pch.list = c(21,22,23,24,25,21,22,23,24,4)
  use.list = c("Pharmaceutical","Herbicide","Pesticide","Organometallic","Antimicrobial","Industrial Chemical","Dye","Antifungal","Fragrance","Other")
  dataset = "heparg2d_toxcast_pfas_pe1_normal"
  celltype ="HepaRG"
  file <- paste0("../output/super_target_boxplot/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,".xlsx")
  print(file)
  heparg = read.xlsx(file)
  mat = heparg
  mat[mat$pod<0.001,"pod"] = 0.001
  mat[mat$burst_pod<0.001,"burst_pod"] = 0.001
  mat = mat[order(mat$pod),]
  plot(mat$pod~mat$burst_pod,type="n",xlab="Burst POD",ylab="POD",xlim=c(0.001,100),ylim=c(0.001,100),log="xy",cex.lab=1.2,cex.axis=1.2,main=celltype,cex=0.5,pch=21,bg="white")
  x0 = 1e-3
  y0 = 1e2
  dy = 2
  for(i in 1:length(col.list)) {
    points(x0,y0,pch=pch.list[i],bg=col.list[i],cex=1)
    text(x0,y0,use.list[i],pos=4)
    y0 = y0/dy
  }
  for(i in 1:nrow(mat)) {
    dtxsid = mat[i,"dtxsid"]
    temp = cann[is.element(cann$dtxsid,dtxsid),]
    use = temp[1,"use_class"]
    if(!is.na(use)) {
      if(is.element(use,use.list)) {
        icol = which.max(use.list==use)
        x = mat[i,"burst_pod"]
        y = mat[i,"pod"]
        delta = abs(rnorm(1,0,0.1))
        if(x<0.001) x = 0.001 * (1+delta)
        if(y<0.001) y = 0.001 * (1+delta)
        points(x,y,pch=pch.list[icol],bg=col.list[icol],cex=1)
      }
    }
  }
  lines(c(0.00001, 100000), c(0.00001,100000))
  lines(c(0.000001,10000),  c(0.00001,100000))
  lines(c(0.0001,  1000000),c(0.00001,100000))

  dataset = "u2os_toxcast_pfas_pe1_normal"
  celltype ="U2OS"
  file <- paste0("../output/super_target_boxplot/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,".xlsx")
  u2os = read.xlsx(file)
  mat = u2os
  plot(mat$pod~mat$burst_pod,type="n",xlab="Burst POD",ylab="POD",xlim=c(0.001,100),ylim=c(0.001,100),log="xy",cex.lab=1.2,cex.axis=1.2,main=celltype,cex=0.5,pch=21,bg="white")
  x0 = 1e-3
  y0 = 1e2
  dy = 2
  for(i in 1:length(col.list)) {
    points(x0,y0,pch=pch.list[i],bg=col.list[i],cex=1)
    text(x0,y0,use.list[i],pos=4)
    y0 = y0/dy
  }
  for(i in 1:nrow(mat)) {
    dtxsid = mat[i,"dtxsid"]
    temp = cann[is.element(cann$dtxsid,dtxsid),]
    use = temp[1,"use_class"]
    if(!is.na(use)) {
      if(is.element(use,use.list)) {
        icol = which.max(use.list==use)
        x = mat[i,"burst_pod"]
        y = mat[i,"pod"]
        delta = abs(rnorm(1,0,0.1))
        if(x<0.001) x = 0.001 * (1+delta)
        if(y<0.001) y = 0.001 * (1+delta)
        points(x,y,pch=pch.list[icol],bg=col.list[icol],cex=1)
      }
    }
  }
  lines(c(0.00001, 100000), c(0.00001,100000))
  lines(c(0.000001,10000),  c(0.00001,100000))
  lines(c(0.0001,  1000000),c(0.00001,100000))

  dataset = "mcf7_ph1_pe1_normal_good_pg"
  celltype ="MCF7"
  file <- paste0("../output/super_target_boxplot/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,".xlsx")
  mcf7 = read.xlsx(file)
  mat = mcf7
  plot(mat$pod~mat$burst_pod,type="n",xlab="Burst POD",ylab="POD",xlim=c(0.001,100),ylim=c(0.001,100),log="xy",cex.lab=1.2,cex.axis=1.2,main=celltype,cex=0.5,pch=21,bg="white")
  x0 = 1e-3
  y0 = 1e2
  dy = 2
  for(i in 1:length(col.list)) {
    points(x0,y0,pch=pch.list[i],bg=col.list[i],cex=1)
    text(x0,y0,use.list[i],pos=4)
    y0 = y0/dy
  }
  for(i in 1:nrow(mat)) {
    dtxsid = mat[i,"dtxsid"]
    temp = cann[is.element(cann$dtxsid,dtxsid),]
    use = temp[1,"use_class"]
    if(!is.na(use)) {
      if(is.element(use,use.list)) {
        icol = which.max(use.list==use)
        x = mat[i,"burst_pod"]
        y = mat[i,"pod"]
        delta = abs(rnorm(1,0,0.1))
        if(x<0.001) x = 0.001 * (1+delta)
        if(y<0.001) y = 0.001 * (1+delta)
        points(x,y,pch=pch.list[icol],bg=col.list[icol],cex=1)
      }
    }
  }
  lines(c(0.00001, 100000), c(0.00001,100000))
  lines(c(0.000001,10000),  c(0.00001,100000))
  lines(c(0.0001,  1000000),c(0.00001,100000))

  dataset = "PFAS_HepaRG"
  celltype ="HepaRG"
  file <- paste0("../output/super_target_boxplot/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,".xlsx")
  print(file)
  heparg_pfas = read.xlsx(file)
  mat = heparg_pfas
  plot(mat$pod~mat$burst_pod,type="n",xlab="Burst POD",ylab="POD",xlim=c(0.001,100),ylim=c(0.001,100),log="xy",cex.lab=1.2,cex.axis=1.2,main=paste("PFAS ",celltype),cex=0.5,pch=21,bg="white")
  x0 = 1e-3
  y0 = 1e2
  dy = 2
  for(i in 1:length(col.list)) {
    points(x0,y0,pch=pch.list[i],bg=col.list[i],cex=1)
    text(x0,y0,use.list[i],pos=4)
    y0 = y0/dy
  }
  for(i in 1:nrow(mat)) {
    dtxsid = mat[i,"dtxsid"]
    temp = cann[is.element(cann$dtxsid,dtxsid),]
    use = temp[1,"use_class"]
    if(!is.na(use)) {
      if(is.element(use,use.list)) {
        icol = which.max(use.list==use)
        x = mat[i,"burst_pod"]
        y = mat[i,"pod"]
        delta = abs(rnorm(1,0,0.1))
        if(x<0.001) x = 0.001 * (1+delta)
        if(y<0.001) y = 0.001 * (1+delta)
        points(x,y,pch=pch.list[icol],bg=col.list[icol],cex=1)
      }
    }
  }
  lines(c(0.00001, 100000), c(0.00001,100000))
  lines(c(0.000001,10000),  c(0.00001,100000))
  lines(c(0.0001,  1000000),c(0.00001,100000))

  dataset = "PFAS_U2OS"
  celltype ="U2OS"
  file <- paste0("../output/super_target_boxplot/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,".xlsx")
  print(file)
  u2os_pfas = read.xlsx(file)
  mat = u2os_pfas
  plot(mat$pod~mat$burst_pod,type="n",xlab="Burst POD",ylab="POD",xlim=c(0.001,100),ylim=c(0.001,100),log="xy",cex.lab=1.2,cex.axis=1.2,main=paste("PFAS ",celltype),cex=0.5,pch=21,bg="white")
  x0 = 1e-3
  y0 = 1e2
  dy = 2
  for(i in 1:length(col.list)) {
    points(x0,y0,pch=pch.list[i],bg=col.list[i],cex=1)
    text(x0,y0,use.list[i],pos=4)
    y0 = y0/dy
  }
  for(i in 1:nrow(mat)) {
    dtxsid = mat[i,"dtxsid"]
    temp = cann[is.element(cann$dtxsid,dtxsid),]
    use = temp[1,"use_class"]
    if(!is.na(use)) {
      if(is.element(use,use.list)) {
        icol = which.max(use.list==use)
        x = mat[i,"burst_pod"]
        y = mat[i,"pod"]
        delta = abs(rnorm(1,0,0.1))
        if(x<0.001) x = 0.001 * (1+delta)
        if(y<0.001) y = 0.001 * (1+delta)
        points(x,y,pch=pch.list[icol],bg=col.list[icol],cex=1)
      }
    }
  }
  lines(c(0.00001, 100000), c(0.00001,100000))
  lines(c(0.000001,10000),  c(0.00001,100000))
  lines(c(0.0001,  1000000),c(0.00001,100000))

  mat1 = heparg
  mat2 = u2os
  dtxsid.list = mat1$dtxsid
  dtxsid.list = dtxsid.list[is.element(dtxsid.list,mat2$dtxsid)]
  plot(c(1,1),type="n",xlab="HepaRG",ylab="U2OS",xlim=c(0.001,100),ylim=c(0.001,100),log="xy",cex.lab=1.2,cex.axis=1.2,main="POD Comparison")
  x0 = 1e-3
  y0 = 1e2
  dy = 2
  for(i in 1:length(col.list)) {
    points(x0,y0,pch=pch.list[i],bg=col.list[i],cex=1)
    text(x0,y0,use.list[i],pos=4)
    y0 = y0/dy
  }
  for(dtxsid in dtxsid.list) {
    x = min(mat1[is.element(mat1$dtxsid,dtxsid),"pod"])
    y = min(mat2[is.element(mat2$dtxsid,dtxsid),"pod"])
    #points(x,y,cex=0.5,pch=21,bg="white")
  }
  for(dtxsid in dtxsid.list) {
    temp = cann[is.element(cann$dtxsid,dtxsid),]
    use = temp[1,"use_class"]
    if(!is.na(use)) {
      if(is.element(use,use.list)) {
        icol = which.max(use.list==use)
        x = min(mat1[is.element(mat1$dtxsid,dtxsid),"pod"])
        y = min(mat2[is.element(mat2$dtxsid,dtxsid),"pod"])
        delta = abs(rnorm(1,0,0.1))
        if(x<0.001) x = 0.001 * (1+delta)
        if(y<0.001) y = 0.001 * (1+delta)
        points(x,y,pch=pch.list[icol],bg=col.list[icol],cex=1)
      }
    }
  }

  lines(c(0.00001, 100000), c(0.00001,100000))
  lines(c(0.000001,10000),  c(0.00001,100000))
  lines(c(0.0001,  1000000),c(0.00001,100000))

  mat1 = heparg_pfas
  mat2 = u2os_pfas
  dtxsid.list = mat1$dtxsid
  dtxsid.list = dtxsid.list[is.element(dtxsid.list,mat2$dtxsid)]
  plot(c(1,1),type="n",xlab="HepaRG PFAS",ylab="U2OS PFAS",xlim=c(0.001,100),ylim=c(0.001,100),log="xy",cex.lab=1.2,cex.axis=1.2,main="POD Comparison")
  x0 = 1e-3
  y0 = 1e2
  dy = 2
  for(i in 1:length(col.list)) {
    points(x0,y0,pch=pch.list[i],bg=col.list[i],cex=1)
    text(x0,y0,use.list[i],pos=4)
    y0 = y0/dy
  }
  for(dtxsid in dtxsid.list) {
    x = min(mat1[is.element(mat1$dtxsid,dtxsid),"pod"])
    y = min(mat2[is.element(mat2$dtxsid,dtxsid),"pod"])
    #points(x,y,cex=0.5,pch=21,bg="white")
  }
  for(dtxsid in dtxsid.list) {
    temp = cann[is.element(cann$dtxsid,dtxsid),]
    use = temp[1,"use_class"]
    if(!is.na(use)) {
      if(is.element(use,use.list)) {
        icol = which.max(use.list==use)
        x = min(mat1[is.element(mat1$dtxsid,dtxsid),"pod"])
        y = min(mat2[is.element(mat2$dtxsid,dtxsid),"pod"])
        delta = abs(rnorm(1,0,0.1))
        if(x<0.001) x = 0.001 * (1+delta)
        if(y<0.001) y = 0.001 * (1+delta)
        points(x,y,pch=pch.list[icol],bg=col.list[icol],cex=1)
      }
    }
  }
  lines(c(0.00001, 100000), c(0.00001,100000))
  lines(c(0.000001,10000),  c(0.00001,100000))
  lines(c(0.0001,  1000000),c(0.00001,100000))
  if(!to.file) browser()

  #par(mfrow=c(3,2),mar=c(4,4,4,4))

  file <- "../toxcast/toxcast_pod.xlsx"
  toxcast = read.xlsx(file)

  mat1 = toxcast
  mat2 = heparg
  dtxsid.list = mat1$dtxsid
  dtxsid.list = dtxsid.list[is.element(dtxsid.list,mat2$dtxsid)]
  plot(c(1,1),type="n",xlab="ToxCast",ylab="HepaRG",xlim=c(0.001,100),ylim=c(0.001,100),log="xy",cex.lab=1.2,cex.axis=1.2,main="POD Comparison")
  x0 = 1e-3
  y0 = 1e2
  dy = 2
  for(i in 1:length(col.list)) {
    points(x0,y0,pch=pch.list[i],bg=col.list[i],cex=1)
    text(x0,y0,use.list[i],pos=4)
    y0 = y0/dy
  }
  for(dtxsid in dtxsid.list) {
    x = min(mat1[is.element(mat1$dtxsid,dtxsid),"pod_uM"])
    y = min(mat2[is.element(mat2$dtxsid,dtxsid),"pod"])
    #points(x,y,cex=0.5,pch=21,bg="white")
  }
  for(dtxsid in dtxsid.list) {
    temp = cann[is.element(cann$dtxsid,dtxsid),]
    use = temp[1,"use_class"]
    if(!is.na(use)) {
      if(is.element(use,use.list)) {
        icol = which.max(use.list==use)
        x = min(mat1[is.element(mat1$dtxsid,dtxsid),"pod_uM"])
        y = min(mat2[is.element(mat2$dtxsid,dtxsid),"pod"])
        delta = abs(rnorm(1,0,0.1))
        if(x<0.001) x = 0.001 * (1+delta)
        if(y<0.001) y = 0.001 * (1+delta)
        points(x,y,pch=pch.list[icol],bg=col.list[icol],cex=1)
      }
    }
  }
  lines(c(0.00001, 100000), c(0.00001,100000))
  lines(c(0.000001,10000),  c(0.00001,100000))
  lines(c(0.0001,  1000000),c(0.00001,100000))

  mat1 = toxcast
  mat2 = mcf7
  dtxsid.list = mat1$dtxsid
  dtxsid.list = dtxsid.list[is.element(dtxsid.list,mat2$dtxsid)]
  plot(c(1,1),type="n",xlab="ToxCast",ylab="MCF7",xlim=c(0.001,100),ylim=c(0.001,100),log="xy",cex.lab=1.2,cex.axis=1.2,main="POD Comparison")
  x0 = 1e-3
  y0 = 1e2
  dy = 2
  for(i in 1:length(col.list)) {
    points(x0,y0,pch=pch.list[i],bg=col.list[i],cex=1)
    text(x0,y0,use.list[i],pos=4)
    y0 = y0/dy
  }
  for(dtxsid in dtxsid.list) {
    x = min(mat1[is.element(mat1$dtxsid,dtxsid),"pod_uM"])
    y = min(mat2[is.element(mat2$dtxsid,dtxsid),"pod"])
    #points(x,y,cex=0.5,pch=21,bg="white")
  }
  for(dtxsid in dtxsid.list) {
    temp = cann[is.element(cann$dtxsid,dtxsid),]
    use = temp[1,"use_class"]
    if(!is.na(use)) {
      if(is.element(use,use.list)) {
        icol = which.max(use.list==use)
        x = min(mat1[is.element(mat1$dtxsid,dtxsid),"pod_uM"])
        y = min(mat2[is.element(mat2$dtxsid,dtxsid),"pod"])
        points(x,y,pch=pch.list[icol],bg=col.list[icol],cex=1)
        delta = abs(rnorm(1,0,0.1))
        if(x<0.001) x = 0.001 * (1+delta)
        if(y<0.001) y = 0.001 * (1+delta)
      }
    }
  }
  lines(c(0.00001, 100000), c(0.00001,100000))
  lines(c(0.000001,10000),  c(0.00001,100000))
  lines(c(0.0001,  1000000),c(0.00001,100000))

  mat1 = toxcast
  mat2 = u2os
  dtxsid.list = mat1$dtxsid
  dtxsid.list = dtxsid.list[is.element(dtxsid.list,mat2$dtxsid)]
  plot(c(1,1),type="n",xlab="ToxCast",ylab="U2OS",xlim=c(0.001,100),ylim=c(0.001,100),log="xy",cex.lab=1.2,cex.axis=1.2,main="POD Comparison")
  x0 = 1e-3
  y0 = 1e2
  dy = 2
  for(i in 1:length(col.list)) {
    points(x0,y0,pch=pch.list[i],bg=col.list[i],cex=1)
    text(x0,y0,use.list[i],pos=4)
    y0 = y0/dy
  }
  for(dtxsid in dtxsid.list) {
    x = min(mat1[is.element(mat1$dtxsid,dtxsid),"pod_uM"])
    y = min(mat2[is.element(mat2$dtxsid,dtxsid),"pod"])
    #points(x,y,cex=0.5,pch=21,bg="white")
  }
  for(dtxsid in dtxsid.list) {
    temp = cann[is.element(cann$dtxsid,dtxsid),]
    use = temp[1,"use_class"]
    if(!is.na(use)) {
      if(is.element(use,use.list)) {
        icol = which.max(use.list==use)
        x = min(mat1[is.element(mat1$dtxsid,dtxsid),"pod_uM"])
        y = min(mat2[is.element(mat2$dtxsid,dtxsid),"pod"])
        delta = abs(rnorm(1,0,0.1))
        if(x<0.001) x = 0.001 * (1+delta)
        if(y<0.001) y = 0.001 * (1+delta)
        points(x,y,pch=pch.list[icol],bg=col.list[icol],cex=1)
      }
    }
  }
  lines(c(0.00001, 100000), c(0.00001,100000))
  lines(c(0.000001,10000),  c(0.00001,100000))
  lines(c(0.0001,  1000000),c(0.00001,100000))

  if(!to.file) browser()
  if(to.file) dev.off()
}

