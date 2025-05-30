#' Pathway Concentration Response Plot
#'
#' Plots a concentration response curve for one sample/signature combination.
#'
#' row is one row of PATHWAY_CR, the signatureConcResp output.
#'
#' @param row Named list containing:
#'   \itemize{
#'     \item conc - conc string separated by |'s
#'     \item resp - response string separated by |'s
#'     \item method - scoring method determines plot bounds
#'     \item proper_name - chemical name for plot title
#'     \item cutoff - noise cutoff
#'     \item bmr - baseline median response; level at which bmd is calculated
#'     \item er - fitted error term for plotting error bars
#'     \item a, tp, b, ga, p, la, q - other model parameters for fit curve
#'     \item fit_method - curve fit method
#'     \item bmd, bmdl, bmdu - bmd, bmd lower bound, and bmd upper bound
#'     \item ac50, acc - curve value at 50% of top, curve value at cutoff
#'     \item top - curve top
#'     \item time, signature, signature_class, signature_size - other identifiers
#'   }
#'   Other elements are ignored.
#' @param plotrange The x-range of the plot as a vector of 2 elements, this can be changed for special cases, but defaults to 0.001 to 100
#' @return No output.
#' @export geneConcRespPlot
#'
#' @importFrom stringr str_split
#' @importFrom grDevices rgb
#' @importFrom graphics lines plot text points axis arrows rect
#' @importFrom stats qt

geneConcRespPlot <- function(row,plotrange=c(0.001,100)) {

  dtxsid <- row[1,"dtxsid"]
  sample_id <- row[1,"sample_id"]

  #every variable in PATHWAY_CR goes into the environment to make it easy
  #to update this function to use new PATHWAY_CR data.
  list2env(row,envir = environment())
  #hard-code plotting points for curves
  logc_plot <- seq(from=-3,to=2,by=0.05)
  conc_plot <- 10**logc_plot

  #reformat conc and resp as vectors
  conc <- as.numeric(str_split(row[1,"conc"],"\\|")[[1]])
  resp <- as.numeric(str_split(row[1,"resp"],"\\|")[[1]])

  #y range is also hard-coded
  ymax <- 3
  ymin <- -3
  #some deprecated code; later will use j =1 and col.list[j] to mean black
  col.list <- c("black","cyan","red")

  #empty plot to start with

  plot(c(1,1),type="n",xlab="conc (uM)",ylab="Score",xlim=plotrange,ylim=c(ymin,ymax),
       log="x",main=paste(proper_name,"\n",gene),cex.main=0.9,cex.lab=1.4,cex.axis=1.4,
       xaxt="n")
  axis(side=1,at=c(0.001,0.01,0.1,1,10,100),cex=1.4)
  delta <- (ymax-ymin)/16
  yval <- ymin

  #cutoffs and gray rectangular noise region
  rect(xleft=0.000001,ybottom=-cutoff,xright=1000,ytop=cutoff,col="lightgray")
  lines(c(0.000001,1000),c(cutoff,cutoff),lwd=1)
  lines(c(0.000001,1000),c(-cutoff,-cutoff),lwd=1)

  #thick line at 0 and bmrs
  lines(c(0.000001,1000),c(0,0),lwd=2)
  lines(c(0.000001,1000),c(bmr,bmr),lwd=1)
  lines(c(0.000001,1000),c(-bmr,-bmr),lwd=1)

  #height for top labels
  yplot <- ymax*0.95
  #equally spaced labelling points pegged to plotrange
  xplot = 10^(seq(log10(plotrange[1]), log10(plotrange[2]), length.out = 8))[-8]

  #Top label headings
  #  text(xplot[1],yplot,"mthd",pos=4)
  text(xplot[1],yplot,"AC50",pos=4)
  text(xplot[2],yplot,"BMD",pos=4)
  text(xplot[3],yplot,"Top",pos=4)
  text(xplot[4],yplot,"T/C",pos=4)
  text(xplot[5],yplot,"Hitcall",pos=4)

  j = 1 #j = 1 is black
  conc <- conc[!is.na(resp)]
  resp <- resp[!is.na(resp)]
  points(resp~conc,pch=19,col=col.list[j]) #plot actual data

  #draw error bars according to er, using 95% confidence for t-dist w/ four df
  if(!is.na(er)){
    arrows(conc, resp+exp(er)*qt(.025,4), conc, resp+exp(er)*qt(.975,4), length = .05, angle = 90, code = 3)
  }

  #get model parameters
  parnames = c("a", "tp", "b", "ga", "p", "la", "q")
  modpars = as.list(row[,parnames])
  modpars= modpars[!sapply(modpars, is.na)]

  #gcalculate and plot model curves
  if(!is.na(fit_method)) {
    if(fit_method == "hill"){
      resp_plot <- do.call("hillfn",list(ps = unlist(modpars), x = conc_plot))
      lines(resp_plot~conc_plot,col=col.list[j])
    } else if(!fit_method %in% c("cnst","none") ){
      resp_plot <- do.call(fit_method,list(ps = unlist(modpars), x = conc_plot))
      lines(resp_plot~conc_plot,col=col.list[j])
    }
  }

  yplot <- yplot-(ymax-ymin)*0.05 #second row for top labels

  #Fill in top labels second row
  text(xplot[1],yplot,format(ac50,digits=2),pos=4, col = "black")
  text(xplot[2],yplot,format(bmd,digits=2),pos=4, col = "black")
  text(xplot[3],yplot,format(top,digits=2),pos=4)
  text(xplot[4],yplot,format(top_over_cutoff,digits=2),pos=4)

  #Bottom left info
  text(xplot[1],0.85*ymin,paste("\nmethod: ",
                                fit_method,"\ncutoff: ",format(cutoff,digits=2),sep=""),pos=4)

  text(1,ymin*0.9,dtxsid,pos=4)
  text(1,ymin,sample_id,pos=4)
  #color hitcall based on whether it's a hit
  color <- "black"
  font <- 1
  if(hitcall==1) {
    color <- "red"
    font <- 2
  }
  text(xplot[5],yplot,format(hitcall, digits = 2),pos=4,col=color,cex=1,font=font)

  #plot green bmd with range
  if(hitcall>0) {
    lines(c(bmd,bmd),c(ymin/2,ymax/2),col="green",lwd=2, lty = isTRUE(bmd<min(conc)) + 1)
    if(is.na(bmdl)) xleft = plotrange[1]/10 else xleft = bmdl
    if(is.na(bmdu)) xright = plotrange[2]*10 else xright = bmdu

    rect(xleft=xleft,ybottom=ymin/2,xright=xright,ytop=ymax/2,col=rgb(0,1,0, alpha = .5), border = NA)
    lines(c(xleft,xleft),c(ymin/2,ymax/2),col="green",lwd=1)
    lines(c(xright,xright),c(ymin/2,ymax/2),col="green",lwd=1)
  }
}
