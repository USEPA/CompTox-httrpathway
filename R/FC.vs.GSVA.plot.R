#--------------------------------------------------------------------------------------
#'
#' Plot the difference between the FC and GSVA analysis at hte signature level
#' @param to.file If TRUE, write plots to a file
#--------------------------------------------------------------------------------------
FC.vs.GSVA.plot <- function(to.file=F,
                            dataset="DMEM_6hr_pilot_normal_pe_0",
                            pathset="PathwaySet_20191107",
                            cutoff= 0.5) {
  printCurrentFunction()
  if(to.file) {
    fname <- paste0("../output/pod_fc_gsva_compare/pod_fc_gsva_compare",dataset,"_",pathset,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,2,2))
  file <- paste0("../output/signature_conc_resp_summary/PATHWAY_CR_",pathset,"_",dataset,"_fc_0.05_conthits.RData")
  print(file)
  load(file)
  fc <- PATHWAY_CR

  file <- paste0("../output/signature_conc_resp_summary/PATHWAY_CR_",pathset,"_",dataset,"_gsva_0.05_conthits.RData")
  print(file)
  load(file)
  gsva <- PATHWAY_CR

  chems <- unique(fc[,c("dtxsid","casrn","name")])
  chems <- chems[order(chems$name),]
  for(i in 1:nrow(chems)) {
    dtxsid <- chems[i,"dtxsid"]
    name <- chems[i,"name"]

    tempfc <- fc[is.element(fc$dtxsid,dtxsid),c("signature","hitcall","bmd")]
    tempgsva <- gsva[is.element(gsva$dtxsid,dtxsid),c("signature","hitcall","bmd")]
    rownames(tempfc) <- tempfc$signature
    rownames(tempgsva) <- tempgsva$signature
    pathlist <- tempfc$signature
    pathlist <- pathlist[is.element(pathlist,tempgsva$signature)]
    tempfc <- tempfc[pathlist,]
    tempgsva <- tempgsva[pathlist,]

    xfc <- tempfc$hitcall
    yfc <- tempfc$bmd
    yfc[xfc<cutoff] <- 1000
    yfc[is.na(yfc)] <- 1000

    xgsva <- tempgsva$hitcall
    ygsva <- tempgsva$bmd
    ygsva[xgsva<cutoff] <- 1000
    ygsva[is.na(ygsva)] <- 1000

    res <- lm(log10(ygsva)~log10(yfc))
    r2 <- summary(res)$adj.r.squared
    rmse <- summary(res)$sigma

    x <- yfc
    y <- ygsva
    x[] <- 0
    y[] <- 0
    x[yfc<1000] <- 1
    y[ygsva<1000] <- 1
    tp <- sum(x*y)
    fp <- sum(x*(1-y))
    fn <- sum((1-x)*y)
    tn <- sum((1-x)*(1-y))

    plot(ygsva~yfc,xlab="POD(FC)",ylab="POD(GSVA)",log="xy",main=name,cex.lab=1.2,cex.axis=1.2,xlim=c(1e-5,1e3),ylim=c(1e-5,1e3))
    text(10**-5, 10**2.5,paste("R2:",format(r2,digits=3)," RMSE:",format(rmse,digits=2)),pos=4)
    text(10**-5, 10**2.0,"tp,fp,fn,tn",pos=4)
    text(10**-5, 10**1.5,paste(tp,fp,fn,tn),pos=4)

    lines(c(1e-10,1e10),c(1e-10,1e10))
    lines(c(1e-11,1e9),c(1e-10,1e10))
    lines(c(1e-9,1e11),c(1e-10,1e10))
    if(!to.file) browser()
  }
  dev.off()
}

