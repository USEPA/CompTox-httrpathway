#----------------------------------------------------------------------------------------
#' Calculate the all-gene cocnentraiton-response and BMDL
#'
#' @param to.file If TRUE, write the plots to a pdf file
#' @param basedir The base directory to find the raw fold change and chemcial data
#' @param dataset The name of the datset to be used
#'
#' @return nothing is returned
#' @export
#----------------------------------------------------------------------------------------
allGeneBMD <- function(to.file=F,basedir="../input/fcdata/",dataset="DMEM_6hr_pilot_normal_00",metric="z",l2fc.limit=1.2) {
  if(to.file) {
    fname <- paste0("../output/signature_conc_resp_summary/AllGeneBMD_ ",dataset," ",metric," ",l2fc.limit,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,2,2))
  ymin <- -1
  ymax <- 1
  file <- "../input/cytotoxicity summary wide allchems.xlsx"
  CYTOTOX <- read.xlsx(file)
  rownames(CYTOTOX) <- CYTOTOX$dtxsid

  aicc <- F
  conthits <- T
  file <- paste0(basedir,"FCMAT2_",dataset,".RData")
  print(file)
  load(file)
  file <- paste0(basedir,"CHEM_DICT_",dataset,".RData",sep="")
  load(file)
  rownames(CHEM_DICT) <- CHEM_DICT[,"sample_key"]

  fcmat <- FCMAT2
  chems <- CHEM_DICT
  fcmat[is.na(fcmat)] <- 0

  chems <- chems[order(chems$name,chems$conc),]
  key.list <- chems$sample_key
  fcmat <- fcmat[key.list,]

  if(metric=="z") {
    cs <- colSums(fcmat)/nrow(fcmat)
    csd <- apply(fcmat,FUN=sd,MARGIN=2)
    z <- sweep(fcmat,2,csd,"/")
    #z <- abs(z)
    meanvec <- apply(z,FUN=mean,MARGIN=1)
    mat <- cbind(chems[,c("dtxsid","casrn","name","conc")],meanvec)
    names(mat)[ncol(mat)] <- "resp"
    ymin <- -1
    ymax <- 1
  }
  if(metric=="zabs") {
    cs <- colSums(fcmat)/nrow(fcmat)
    csd <- apply(fcmat,FUN=sd,MARGIN=2)
    z <- sweep(fcmat,2,csd,"/")
    z <- abs(z)
    meanvec <- apply(z,FUN=mean,MARGIN=1)
    mat <- cbind(chems[,c("dtxsid","casrn","name","conc")],meanvec)
    names(mat)[ncol(mat)] <- "resp"
    ymin <- 0
    ymax <- 2
  }
  if(metric=="l2fc") {
    cs <- colSums(fcmat)/nrow(fcmat)
    meanvec <- apply(fcmat,FUN=mean,MARGIN=1)
    mat <- cbind(chems[,c("dtxsid","casrn","name","conc")],meanvec)
    names(mat)[ncol(mat)] <- "resp"
    ymin <- -0.2
    ymax <- 0.2
  }
  if(metric=="deg") {
    temp <- abs(fcmat)
    temp[temp<l2fc.limit] <- 0
    temp[temp>0] <- 1
    meanvec <- rowSums(temp)/ncol(temp)
    mat <- cbind(chems[,c("dtxsid","casrn","name","conc")],meanvec)
    names(mat)[ncol(mat)] <- "resp"
    ymin <- 0
    ymax <- 0.5
  }

  dtxsid.list <- unique(mat$dtxsid)
  name.list <- c("dtxsid","casrn","name","bmd","bmdu","bmdl","top","hitcall")
  result <- as.data.frame(matrix(nrow=length(dtxsid.list),ncol=length(name.list)))
  names(result) <- name.list
  counter <- 0
  temp <- sort(mat$resp)
  counter <- floor(0.66*length(temp))
  onesd <- sd(mat$resp)
  for(dtxsid in dtxsid.list) {
    #dtxsid <- "DTXSID6024882"
    temp <- mat[is.element(mat$dtxsid,dtxsid),]
    casrn <- temp[1,"casrn"]
    name <- temp[1,"name"]
    counter <- counter+1
    result[counter,"dtxsid"] <- dtxsid
    result[counter,"casrn"] <- casrn
    result[counter,"name"] <- name
    cat(name,"\n")
    resp <- temp$resp
    conc <- temp$conc
    logc <- log10(conc)
    cutoff <- onesd
    params <- httrFit(conc=conc, resp=resp, cutoff=cutoff, force.fit = T, bidirectional = TRUE, verbose = F, do.plot = F,
                      fitmodels = c("cnst", "hill", "poly1", "poly2", "pow", "exp2", "exp3", "exp4", "exp5","gnls"))

    ##################################################

    #initialize parameters to NA
    a = b = tp = p = q = ga = la = er = top = ac50 = ac50_loss = ac5 = ac10 = ac20 = acc = ac1sd = bmd = NA_real_
    bmdl = bmdu = caikwt = mll = NA_real_

    #get aics and degrees of freedom
    aics = sapply(params$modelnames, function(x){params[[x]][["aic"]]})
    dfs = sapply(params$modelnames, function(x){length(params[[x]][["pars"]])})
    aics = aics[!is.na(aics)]

    if(sum(!is.na(aics)) == 0){
      #if all fits failed, use none for method
      fit_method = "none"
      rmse = NA_real_
    } else {
      #use nested chisq to choose between poly1 and poly2, remove poly2 if it fails.
      #pvalue hardcoded to .05
      aics = nestselect(aics, "poly1", "poly2", dfdiff = 1, pval = .05)
      dfs = dfs[names(dfs) %in% names(aics)]

      #it's useful to keep original aics so we can extract loglikelihoods for nested models (above) and mll (below)
      if(aicc)   saics = aics + 2*dfs*(dfs+1)/(length(resp)-dfs-1)   else saics = aics

      if(conthits) {
        #get aikaike weight of winner (vs constant) for cont hitcalls
        #never choose constant as winner for cont hitcalls
        nocnstaics = saics[names(saics) != "cnst"]
        fit_method = names(nocnstaics)[which.min(nocnstaics)]
        caikwt = exp(-saics["cnst"]/2)/(exp(-saics["cnst"]/2) + exp(-saics[fit_method]/2))
        caikwt = 1 + exp(-(saics[fit_method]-saics["cnst"])/2)
        caikwt = 1/ caikwt
        #cat("caikwt",caikwt,"\n")
        if(is.nan(caikwt)) caikwt <- 1
      } else  {
        fit_method = names(saics)[which.min(saics)]
      }

      fitout = params[[fit_method]]
      rmse = fitout$rme
      modpars = fitout[fitout$pars]
      list2env(modpars, envir = environment()) #put model parameters in environment
    }

    #first deal with parameter output
    if(fit_method %in% c("poly1", "poly2", "pow", "exp2", "exp3")){
      #methods that grow without bound: top defined as model value at max conc
      top = fitout$modl[which.max(abs(fitout$modl))] #top is taken to be highest model value
      ac50 = acy(.5*top, modpars, type = fit_method)
    } else if(fit_method %in% c("hill", "exp4", "exp5")){
      #methods with a theoretical top/ac50
      top = tp
      ac50 = ga
    } else if(fit_method == "gnls"){
      #gnls methods; use calculated top/ac50, etc.
      top =  acy(0, modpars, type = fit_method, returntop = T)
      ac50 = acy(.5*top, modpars, type = fit_method)
      ac50_loss = acy(.5*top, modpars, type = fit_method, getloss = T)
    }
    n_gt_cutoff = sum(abs(resp)>cutoff)

    #compute discrete or continuous hitcalls
    if(fit_method == "none") {
      hitcall = 0
    } else if(conthits){
      mll = length(modpars) - aics[[fit_method]]/2
      hitcall = hitcontinner(conc, resp, top, cutoff, er, ps = unlist(modpars), fit_method,
                             caikwt = caikwt, mll = mll)
    } else {
      hitcall = hitloginner(conc, resp, top, cutoff, ac50)
    }
    bmr = onesd*1.349 #magic bmr is hard-coded
    if(hitcall > 0){
      #fill ac's; can put after hit logic
      ac5 = acy(.05*top, modpars, type = fit_method) #note: cnst model automatically returns NAs
      ac10 = acy(.1*top, modpars, type = fit_method)
      ac20 = acy(.2*top, modpars, type = fit_method)
      acc = acy(sign(top)*cutoff, modpars, type = fit_method)
      ac1sd = acy(sign(top)*onesd, modpars, type = fit_method)
      bmd = acy(sign(top)*bmr, modpars, type = fit_method)

      #get bmdl and bmdu
      bmdl = bmdbounds(fit_method, bmr = sign(top)*bmr, pars = unlist(modpars), conc, resp, onesidedp = .05,
                       bmd = bmd, which.bound = "lower")
      bmdu = bmdbounds(fit_method, bmr = sign(top)*bmr, pars = unlist(modpars), conc, resp, onesidedp = .05,
                       bmd = bmd, which.bound = "upper")
    }
    top_over_cutoff <- abs(top)/cutoff
    conc.plot <- paste(conc,collapse="|")
    resp.plot <- paste(resp,collapse="|")

    #PATHWAY_CR contains the specified columns and any identifying, unused columns
    #that were in signaturescoremat/
    #identifiers = row[!names(row) %in% c("conc", "resp", "bmed", "onesd", "cutoff")]
    name.list <- c("onesd", "cutoff",
                   "n_gt_cutoff","cutoff", "fit_method",
                   "top_over_cutoff", "rmse", "a", "b", "tp", "p", "q", "ga", "la", "er", "bmr", "bmdl", "bmdu", "caikwt",
                   "mll","hitcall", "ac50","ac50_loss","top", "ac5","ac10","ac20", "acc", "ac1sd", "bmd", "conc.plot", "resp.plot")
    row = as.data.frame(c(mget(name.list)), stringsAsFactors = F)

    ##################################################

    bmd <- log10(bmd)
    bmdl <- log10(bmdl)
    bmdu <- log10(bmdu)

    #cat("bmd: ",row$bmd,"\n")
    plot(resp~logc,main=paste(name),cex.lab=1.2,cex.axis=1.2,xlab="log(conc uM)",ylab=paste("Fraction L2fc>",l2fc.limit),
         ylim=c(ymin,ymax),xlim=c(-2,2),pch=21,cex=1,bg="black")

    luc <- CYTOTOX[dtxsid,"LUC"]
    bla <- CYTOTOX[dtxsid,"BLA"]
    srb <- CYTOTOX[dtxsid,"SRB"]
    other <- CYTOTOX[dtxsid,"LUC"]
    delta <- (ymax-ymin)/8
    yval <- ymin
    if(!is.na(bla)) {lines(c(bla,bla),c(yval,yval+delta),lwd=4,col="red"); yval <- yval+delta/2}
    if(!is.na(luc)) {lines(c(luc,luc),c(yval,yval+delta),lwd=4,col="orange"); yval <- yval+delta/2}
    if(!is.na(srb)) {lines(c(srb,srb),c(yval,yval+delta),lwd=4,col="cyan"); yval <- yval+delta/2}
    if(!is.na(other)) {lines(c(other,other),c(yval,yval+delta),lwd=4,col="blue"); yval <- yval+delta/2}


     #get model parameters
    parnames = c("a", "tp", "b", "ga", "p", "la", "q")
    modpars = as.list(row[,parnames])
    modpars= modpars[!sapply(modpars, is.na)]
    #hard-code plotting points for curves
    logc_plot <- seq(from=-3,to=2,by=0.05)
    conc_plot <- 10**logc_plot

    #calculate and plot model curves
    if(fit_method == "hill"){
      resp_plot <- do.call("hillfn",list(ps = unlist(modpars), x = conc_plot))
      lines(resp_plot~logc_plot)
    } else if(!fit_method %in% c("cnst","none") ){
      resp_plot <- do.call(fit_method,list(ps = unlist(modpars), x = conc_plot))
      lines(resp_plot~logc_plot)
    }

    plotrange <- c(-3,2)
    lines(c(bmd,bmd),c(ymin/2,ymax/2),col="green",lwd=2, lty = isTRUE(bmd<min(conc)) + 1)
    if(is.na(bmdl)) xleft = plotrange[1]/10 else xleft = bmdl
    if(is.na(bmdu)) xright = plotrange[2]*10 else xright = bmdu

    if(hitcall>0.2) {
      rect(xleft=xleft,ybottom=ymin/2,xright=xright,ytop=ymax/2,col=rgb(0,1,0, alpha = .5), border = NA)
      lines(c(xleft,xleft),c(ymin/2,ymax/2),col="green",lwd=1)
      lines(c(xright,xright),c(ymin/2,ymax/2),col="green",lwd=1)
    }

    bmdl.text <- ">100 uM"
    if(!is.na(bmdl)) bmdl.text <- paste(format(10**(bmdl),digits=2),"uM")
    text(-2,0.9*ymax,paste("Hitcall:",format(hitcall,digits=2)),pos=4)
    text(-2,0.8*ymax,paste("BMDL:",bmdl.text),pos=4)
    lines(c(-3,3),c(0,0))
    #limit1 <- cytotox.old[dtxsid,"cytotox_median_log"]
    #limit2 <- log10(cytotox.old[dtxsid,"cytotox_lower_bound_um"])
    #lines(c(limit1,limit1),c(-50,50),col="blue",lwd=2)
    #lines(c(limit2,limit2),c(-50,50),col="red",lwd=2)

    result[counter,"fit_method"] <- fit_method
    result[counter,"bmd"] <- bmd
    result[counter,"bmdu"] <- bmdu
    result[counter,"bmdl"] <- bmdl
    result[counter,"top"] <- top
    result[counter,"hitcall"] <- hitcall
    print(result[counter,])
    if(!to.file) browser()
  }
  file <- paste0("../output/signature_conc_resp_summary/AllGeneBMD_ ",dataset," ",metric," ",l2fc.limit,".xlsx")
  write.xlsx(result,file)
  if(to.file) dev.off()
}
