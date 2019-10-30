#' BMD Accumulation Plot With Nulls
#'
#' Creates pathway BMD accumulation plot vs. null and computes accumulation BMD.
#'
#' Nullset and dataset should already have been run through pathwayConcResp using
#' the given pathset, method, and oldpval. Only generates proof of concept plots
#' for accumulation BMDs. There is not currently a method to extract the
#' accumulation BMDs directly.
#'
#' @param pathset Name of pathset.
#' @param dataset Name of dataset.
#' @param method Name of pathway scoring method.
#' @param nullset Name of null dataset.
#' @param newpval P-value for cutoff to be used in plot.
#' @param oldpval P-value that nullset and dataset were originally run with.
#' @param to.file If to.file = T, plots to file.
#' @param usecont Set usecont = T for continuous hitcalls, usecont = F for discrete.
#'   Should probably match the original hitcall type use in CR modeling.
#' @param nametag Set additional name descriptor that was attached to CR modeling,
#'   if applicable.
#' @param mc.cores Number of cores to use when altering continuous hitcalls; has
#'   no effect if usecont = F or newpval = oldpval.
#'
#' @import grDevices
#'
#' @return No output.
#' @export
pathwayAccumNullPlot = function(pathset="bhrr", dataset="arer", method = "fc", nullset = "arer_RAND125",
                                newpval = .05, oldpval = .05, to.file = T, usecont = T, nametag = NULL, mc.cores = 1){

  #get pathway_Crs from regular data
  if(!is.null(nametag)) nametag = paste0("_",nametag)
  file <- paste0("../output/pathway_conc_resp_summary/PATHWAY_CR_" ,pathset,"_",dataset,"_",method , "_", oldpval,
                 nametag,".RData")
  load(file)

  #remove missing bmd's and get cutoff for newpval
  PATHWAY_CR = PATHWAY_CR[!is.na(PATHWAY_CR$bmd10),]
  pvalkey = getpvalcutoff(pathset, nullset = nullset, method = method, pvals = newpval)
  newcutoff = pvalkey$cutoff[match(PATHWAY_CR$pathway, pvalkey$pathway)]

  #recalc hitcalls if newpval does not equal oldpval
  if(usecont){
    if(oldpval != newpval) PATHWAY_CR$hitcall = hitcont(PATHWAY_CR, newcutoff = newcutoff, mc.cores = mc.cores)
    hitmat = PATHWAY_CR
  } else {
    if(oldpval != newpval) PATHWAY_CR$hitcall = hitlogic(PATHWAY_CR, newcutoff = newcutoff)
    hitmat = PATHWAY_CR[PATHWAY_CR$hitcall == 1,]
  }

  #get bottom, top, chems, order output by chemical name
  bottom = floor(min(log10(hitmat$bmd10 ) ) )
  top = 2
  hitmat = hitmat[order(hitmat$name),]
  sids = unique(hitmat$sample_id)

  #plotting points at every unique bmd up to two significant digits.
  pts = sort(unique(signif(hitmat$bmd10,2)))

  #set colors
  cols = c("black", "blue", "red")

  #get null info
  file <- paste0("../output/pathway_conc_resp_summary/PATHWAY_CR_",pathset,"_",nullset,"_",method ,"_",oldpval, nametag,
                 ".RData")
  load(file)
  PATHWAY_CR = PATHWAY_CR[!is.na(PATHWAY_CR$bmd10),]

  #recalc null hitcalls, if necessary
  pvalkey = getpvalcutoff(pathset, nullset = nullset, method = method, pvals = newpval)
  newcutoff = pvalkey$cutoff[match(PATHWAY_CR$pathway, pvalkey$pathway)]
  if(usecont){
    if(oldpval != newpval) PATHWAY_CR$hitcall = hitcont(PATHWAY_CR, newcutoff = newcutoff, mc.cores = mc.cores)
    nullmat = PATHWAY_CR
  } else {
    if(oldpval != newpval) PATHWAY_CR$hitcall = hitlogic(PATHWAY_CR, newcutoff = newcutoff)
    nullmat = PATHWAY_CR[PATHWAY_CR$hitcall == 1,]
  }

  #this line is probably deprecated, not sure.
  PATHWAY_CR$method = rep(method, nrow(PATHWAY_CR))

  nullnames = unique(nullmat$name)

  #smooth ecdf for each null chem
  outmat = sapply(nullnames, function(chem){
    mymat = nullmat[nullmat$name == chem, ]
    if(nrow(mymat) > 0) return(smoothecdf(pts, mymat)) else return(rep(0, length(pts)))
  }, simplify = "array")

  #take mean, upper 95%, lower 5% of null chems at each pt
  quant =  apply(outmat,1, mean)
  quantu =  apply(outmat,1, quantile, probs = .05)
  quantl =  apply(outmat,1, quantile, probs = .95)

  maxnull = max(quantl)

  #open pdf
  if(to.file){
    dir.create("../output/accumplots/", showWarnings = F)
    file = paste0("../output/accumplots/NULLACCUMPLOT_",pathset,"_",dataset,"_", newpval, nametag, ".pdf")
    pdf(file= file,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    par(mfrow=c(3,2),mar=c(4,4,2,2))
  }

  #get chem dict to show min conc for each chem
  load(paste0("../input/fcdata/CHEM_DICT_",dataset,".RData"))

  #plotting loops, outer is chems
  out = lapply(as.list(sids), function(x){

    mymat = hitmat[hitmat$sample_id == x,]

    #get max number of hits for a chem
    maxhits = max(c(sum(mymat$hitcall), maxnull) )

    #empty plot
    plot(0, type = 'n', xlim = c(bottom, top), ylim = c(0,maxhits),
         ylab = "Pathway Hits (Smoothed)", xlab = "Concentration", main = mymat$name[1], xaxt = "n")
    axis(side = 1, at = c(-6, -5, -4, -3,-2,-1,0,1,2), labels = c(1e-6, 1e-5, 1e-4,.001, .01, .1, 1, 10, 100))

    #draw min conc line
    minconc = min(CHEM_DICT$conc[CHEM_DICT$name == x])
    lines(x = rep(log10(minconc), 2),y= c(0, maxhits), lty = 2)

    if(nrow(mymat) > 0){
      #get and plot ecdf
      ys = smoothecdf(pts, mymat)
      points(log10(pts), ys, type = "l")

      #plot null mean, find and plot accum bmd
      points(log10(pts) , quant, type = "l", col = "red")
      loc = pts[ min(which(ys > (quant + .1*(max(quant))) )) ]
      lines(x= rep(log10(loc),2), y= c(0, maxhits), lwd = 2, col = "blue" )

      #find and plot accum bmdu, bmdl, and add shading between
      locl = pts[ min(which(ys > (quantl + .1*(max(quantl))) )) ]
      if(is.na(locl)) locl = max(pts)
      locu = pts[ min(which(ys > (quantu + .1*(max(quantu))) )) ]
      if(is.na(locu)) locu = max(pts)
      rect(log10(locl), 0, log10(locu), maxhits, col = rgb(0,0,1,.2), border = NA)
      lines(x= rep(log10(locl),2), y= c(0, maxhits), lwd = 1, col = "blue" )
      lines(x= rep(log10(locu),2), y= c(0, maxhits), lwd = 1, col = "blue" )

      #plot null upper/lower dotted lines
      points(log10(pts) , quantl, type = "l", col = "red", lty = 2)
      points(log10(pts) , quantu, type = "l", col = "red", lty = 2)

    }

    legend("topleft", legend = c("BMD Smooth ECDF", "Null Mean", "Null 90% CI",
                                 "Overall BMD", "Overall BMD 90% CI", "Min Concentration"),
           col = c("black", "red", "red", "blue", NA, "black"), cex = .8, lty = c(1,1,2,1,NA, 2),
           fill = c(NA,NA,NA,NA,rgb(0,0,1,.2), NA), border = c(NA,NA,NA,NA,NA, NA))
  })
  if(to.file) dev.off()

}

#' Smooth ECDF
#'
#' Converts a data frame containing bmd10, bmdl, bmu, to a smooth ecdf.
#'
#' Models each bmd as a gaussian with mean bmd10 uses bmdl (bmdu if bmdl is na)
#' to compute sd. Each gaussian is scaled by the hitcall.
#'
#' @param x ECDF plotting location x-values.
#' @param mymat Dataframe containing bmd10, bmdu, bmdl, and hitcall columns.
#' @param verbose verbose = F suppresses both bounds NA warning.
#' @param bmdrange Maximum expected BMD range. The farthest value from the bmd10
#'  is used to compute standard deviation of gaussian when both bounds are missing.
#' @return Outputs a vector corresponding to the locations in x.
#' @export
#'
#' @examples
#' x = 10^(-50:50/30)
#' mymat =data.frame(list(bmd10 = c(.1,1,10), bmdl = c(.05,NA,NA),
#'   bmdu = c(.6,1.5,NA), hitcall = c(1,1,1)))
#' out = smoothecdf(x, mymat)
#' plot(log10(x),out, type = "l")
#' mymat$hitcall = c(1,.5,0)
#' out2 = smoothecdf(x, mymat)
#' plot(log10(x),out2, type = "l")
smoothecdf = function(x, mymat, verbose = F, bmdrange = c(.001,100)){

  output = rep(0,length(x))
  for(i in 1:nrow(mymat)){
    row = mymat[i,]
    mymean = row$bmd10
    if(is.na(mymean)) next #skip if bmd is NA
    #try to use bmdl, then try bmdu, then use bmdrange to get gaussian sd
    if(!is.na(row$bmdl)){
      mysd = (row$bmdl - mymean)/qnorm(.05)
    } else if(!is.na(row$bmdu)){
      mysd = (row$bmdu - mymean)/qnorm(.95)
    } else {
      if(verbose) warning("Both bounds NA, using max range")
      if(mymean < 10^(mean(log10(bmdrange))) ){
        #use upper range as bmdu if bmd less than half of log space range, else use lower range as bmdl
        mysd = (bmdrange[2] - mymean)/qnorm(.95)
      } else {
        mysd = (bmdrange[1] - mymean)/qnorm(.05)
      }
    }

    #add gaussian cdf scaled by hitcall
    output = output + pnorm(x,mean = mymean, sd = mysd)*row$hitcall

  }
  return(output)

}


