#--------------------------------------------------------------------------------------

#'
#' u2os_pilot_pe1_normal_null_full
#' u2os_pilot_pe1_normal_null_pilot
#' u2os_pilot_pe1_normal_null_pilot_lowconc
#' u2os_pilot_pe1_normal_null_pilot_lowconc_lowchem
#'

#'
#'
#--------------------------------------------------------------------------------------
u2os_pilot_null_compare <- function(to.file=F,
                                    sigset="screen_large",
                                    method="fc") {

  printCurrentFunction()


  dataset = "u2os_pilot_pe1_normal_null_full"
  file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
  print(file)
  load(file=file)
  null_full = SIGNATURE_CR
  rownames(null_full) = paste0(null_full$sample_id,"_",null_full$signature)

  dataset = "u2os_pilot_pe1_normal_null_pilot"
  file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
  load(file=file)
  null_pilot = SIGNATURE_CR
  rownames(null_pilot) = paste0(null_pilot$sample_id,"_",null_pilot$signature)

  dataset = "u2os_pilot_pe1_normal_null_pilot_lowconc"
  file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
  load(file=file)
  null_pilot_lowconc = SIGNATURE_CR
  rownames(null_pilot_lowconc) = paste0(null_pilot_lowconc$sample_id,"_",null_pilot_lowconc$signature)

  dataset = "u2os_pilot_pe1_normal_null_pilot_lowconc_lowchem"
  file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
  load(file=file)
  null_pilot_lowconc_lowchem = SIGNATURE_CR
  rownames(null_pilot_lowconc_lowchem) = paste0(null_pilot_lowconc_lowchem$sample_id,"_",null_pilot_lowconc_lowchem$signature)

  null_pilot = null_pilot[rownames(null_full),]
  null_pilot_lowconc = null_pilot_lowconc[rownames(null_full),]
  null_pilot_lowconc_lowchem = null_pilot_lowconc_lowchem[rownames(null_full),]


  ###############################################################################
  if(to.file) {
    fname <- paste0("../output/u2os_pilot/u2os_pilot_null_compare.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }

  ###############################################################################
  par(mfrow=c(3,2),mar=c(4,4,2,2))

  nx = null_full
  ny = null_pilot
  ny = ny[nx$hitcall>=0.9,]
  nx = nx[nx$hitcall>=0.9,]
  x = nx$bmd
  y = ny$bmd
  plot(y~x,cex.lab=1.2,cex.axis=1.2,main="BMD",pch=".",log="xy",
       xlim=c(0.001,1000),ylim=c(0.001,1000),xlab="NULL Full",ylab="NULL Pilot All")
  lines(c(0,1),c(0,1))
  lx = log10(x)
  ly = log10(y)
  res = lm(ly~lx)
  ss=summary(res)
  text(1,0.01,paste("R2=",format(ss$r.squared,digits=2)),pos=4,cex=1.5)
  text(1,0.06,paste("slope=",format(ss$coefficients[2,1],digits=2)),pos=4,cex=1.5)
  lines(c(1e-10,1e10),c(1e-10,1e10))

  nx = null_full
  ny = null_pilot_lowconc
  ny = ny[nx$hitcall>=0.9,]
  nx = nx[nx$hitcall>=0.9,]
  x = nx$bmd
  y = ny$bmd
  plot(y~x,cex.lab=1.2,cex.axis=1.2,main="BMD",pch=".",log="xy",
       xlim=c(0.001,1000),ylim=c(0.001,1000),xlab="NULL Full",ylab="NULL Pilot LowConc")
  lines(c(0,1),c(0,1))
  lx = log10(x)
  ly = log10(y)
  res = lm(ly~lx)
  ss=summary(res)
  text(1,0.01,paste("R2=",format(ss$r.squared,digits=2)),pos=4,cex=1.5)
  text(1,0.06,paste("slope=",format(ss$coefficients[2,1],digits=2)),pos=4,cex=1.5)
  lines(c(1e-10,1e10),c(1e-10,1e10))

  nx = null_full
  ny = null_pilot_lowconc_lowchem
  ny = ny[nx$hitcall>=0.9,]
  nx = nx[nx$hitcall>=0.9,]
  x = nx$bmd
  y = ny$bmd
  plot(y~x,cex.lab=1.2,cex.axis=1.2,main="BMD",pch=".",log="xy",
       xlim=c(0.001,1000),ylim=c(0.001,1000),xlab="NULL Full",ylab="NULL Pilot LowConc LowChem")
  lines(c(0,1),c(0,1))
  lx = log10(x)
  ly = log10(y)
  res = lm(ly~lx)
  ss=summary(res)
  text(1,0.01,paste("R2=",format(ss$r.squared,digits=2)),pos=4,cex=1.5)
  text(1,0.06,paste("slope=",format(ss$coefficients[2,1],digits=2)),pos=4,cex=1.5)
  lines(c(1e-10,1e10),c(1e-10,1e10))
  if(!to.file) browser()

  ###############################################################################
  par(mfrow=c(3,2),mar=c(4,4,2,2))

  x = null_full$hitcall
  y = null_pilot$hitcall
  y = y[x>=0.9]
  x = x[x>=0.9]
  plot(y~x,cex.lab=1.2,cex.axis=1.2,main="Hitcalls",pch=".",
       xlim=c(0,1),ylim=c(0,1),xlab="NULL Full",ylab="NULL Pilot All")
  lines(c(0,1),c(0,1))
  res = lm(y~x)
  ss=summary(res)
  text(0.1,0.8,paste("R2=",format(ss$r.squared,digits=2)),pos=4,cex=1.5)
  text(0.1,0.6,paste("slope=",format(ss$coefficients[2,1],digits=2)),pos=4,cex=1.5)
  plot(density(y),cex.lab=1.2,cex.axis=1.2,main="Hitcalls",xlab="NULL Pilot All")
  n1 = length(y)
  n2 = length(y[y>=0.9])
  ratio = n2/n1
  text(0.1,5,paste("f(>=0.9)=",format(ratio,digits=2)),pos=4,cex=1.5)

  x = null_full$hitcall
  y = null_pilot_lowconc$hitcall
  y = y[x>=0.9]
  x = x[x>=0.9]
  plot(y~x,cex.lab=1.2,cex.axis=1.2,main="Hitcalls",pch=".",
       xlim=c(0,1),ylim=c(0,1),xlab="NULL Full",ylab="NULL Pilot LowConc")
  lines(c(0,1),c(0,1))
  res = lm(y~x)
  ss=summary(res)
  text(0.1,0.8,paste("R2=",format(ss$r.squared,digits=2)),pos=4,cex=1.5)
  text(0.1,0.6,paste("slope=",format(ss$coefficients[2,1],digits=2)),pos=4,cex=1.5)
  plot(density(y),cex.lab=1.2,cex.axis=1.2,main="Hitcalls",xlab="NULL Pilot LowConc")
  n1 = length(y)
  n2 = length(y[y>=0.9])
  ratio = n2/n1
  text(0.1,10,paste("f(>=0.9)=",format(ratio,digits=2)),pos=4,cex=1.5)

  x = null_full$hitcall
  y = null_pilot_lowconc_lowchem$hitcall
  y = y[x>=0.9]
  x = x[x>=0.9]
  plot(y~x,cex.lab=1.2,cex.axis=1.2,main="Hitcalls",pch=".",
       xlim=c(0,1),ylim=c(0,1),xlab="NULL Full",ylab="NULL Pilot LowConc LowChem")
  lines(c(0,1),c(0,1))
  res = lm(y~x)
  ss=summary(res)
  text(0.1,0.8,paste("R2=",format(ss$r.squared,digits=2)),pos=4,cex=1.5)
  text(0.1,0.6,paste("slope=",format(ss$coefficients[2,1],digits=2)),pos=4,cex=1.5)
  plot(density(y),cex.lab=1.2,cex.axis=1.2,main="Hitcalls",xlab="NULL Pilot LowConc LowChem")
  n1 = length(y)
  n2 = length(y[y>=0.9])
  ratio = n2/n1
  text(0.1,20,paste("f(>=0.9)=",format(ratio,digits=2)),pos=4,cex=1.5)
  if(!to.file) browser()

  ###############################################################################
  par(mfrow=c(3,2),mar=c(4,4,2,2))

  plot(null_pilot$cutoff~null_full$cutoff,cex.lab=1.2,cex.axis=1.2,main="Cutoffs",pch=".",
       xlim=c(0,0.2),ylim=c(0,0.2),xlab="NULL Full",ylab="NULL Pilot All")
  lines(c(0,1),c(0,1))
  res = lm(null_pilot$cutoff~null_full$cutoff)
  ss=summary(res)
  text(0.12,0.05,paste("R2=",format(ss$r.squared,digits=2)),pos=4,cex=1.5)
  text(0.12,0.02,paste("slope=",format(ss$coefficients[2,1],digits=2)),pos=4,cex=1.5)

  plot(null_pilot_lowconc$cutoff~null_full$cutoff,cex.lab=1.2,cex.axis=1.2,main="Cutoffs",pch=".",
       xlim=c(0,0.2),ylim=c(0,0.2),xlab="NULL Full",ylab="NULL Pilot LowConc")
  lines(c(0,1),c(0,1))
  res = lm(null_pilot_lowconc$cutoff~null_full$cutoff)
  ss=summary(res)
  text(0.12,0.05,paste("R2=",format(ss$r.squared,digits=2)),pos=4,cex=1.5)
  text(0.12,0.02,paste("slope=",format(ss$coefficients[2,1],digits=2)),pos=4,cex=1.5)

  plot(null_pilot_lowconc_lowchem$cutoff~null_full$cutoff,cex.lab=1.2,cex.axis=1.2,main="Cutoffs",pch=".",
       xlim=c(0,0.2),ylim=c(0,0.2),xlab="NULL Full",ylab="NULL Pilot LowConc LowChem")
  lines(c(0,1),c(0,1))
  res = lm(null_pilot_lowconc_lowchem$cutoff~null_full$cutoff)
  ss=summary(res)
  text(0.12,0.05,paste("R2=",format(ss$r.squared,digits=2)),pos=4,cex=1.5)
  text(0.12,0.02,paste("slope=",format(ss$coefficients[2,1],digits=2)),pos=4,cex=1.5)

  if(!to.file) browser()
  else dev.off()
}
