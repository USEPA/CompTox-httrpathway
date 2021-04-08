#' Plot Outer
#'
#' Calls signatureConcResp plotting function.
#'
#' Calls signatureConcResp plotting function for one chemical and every signature.
#' Saves a single pdf to disk for the given chemical containing every signature
#' CR plot.
#'
#' @param proper_name Chemical name to be used in file name.
#' @param GENE_CR Dataframe output of geneConcResp_pval.
#' @param foldname Folder name for output file.
#' @param plotrange The x-range of the plot as a vector of 2 elements, this can be changed for special cases, but defaults to 0.001 to 100
#' @import grDevices
#'
#' @return No output.
#' @export
plotouterGene = function(proper_name, GENE_CR, foldname,plotrange=c(0.001,100)){
  #open pdf for plots
  #printCurrentFunction(proper_name)
  temp <- GENE_CR[GENE_CR$proper_name == proper_name,]
  sid.list <- unique(temp$sample_id)
  for(sid in sid.list) {
    fname <- paste0(foldname,"/conc_resp_",proper_name," ",sid,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    par(mfrow=c(3,2),mar=c(4,4,2,2))
    par(mfrow=c(3,2),mar=c(5,5,4,2))

    #narrow down to given chemical / and sample_id
    subframe = temp[is.element(temp$sample_id,sid),]

    subframe$auc <- abs(subframe$top) * (3-log10(subframe$bmd))
    subframe <- subframe[order(subframe$auc,decreasing=T),]

    #subframe = subframe[order(-subframe$hitcall,-subframe$top_over_cutoff, subframe$bmd),] #order by potency (optional)

    #cycle through signatures (rows) and run signatureConcRespPlot
    for(i in 1:nrow(subframe)){
      geneConcRespPlot(subframe[i,],plotrange=plotrange)
    }
    dev.off()
  }
  graphics.off()
}


