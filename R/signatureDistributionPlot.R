#' Pathway Distribution Plot
#'
#' Plots null and actual pdfs for given signature and cutoffs.
#'
#' This function requires that a signaturescoremat file has already been generated
#' for the given sigset/dataset/method using signatureScore. There should also
#' be signaturescoremat file for the nullset if comparetype = "Null". This function
#' has also been used to get crossing-based cutoffs, but that feature has been
#' deprecated.
#'
#' @param sigset Name of signature set.
#' @param dataset Name of data set.
#' @param method Name of signature scoring method.
#' @param nullset Name of null data set.
#' @param perc 1-p-value for pvalue cutoff.
#' @param fdr False discovery rate for FDR cutoff.
#' @param comparetype Type of noise to use: "Null" for null data scores,
#'   "Low Conc" for lowest concentrations.
#' @param samplepaths Vector of sample signature names to plot.
#' @param to.file If to.file = T, write plot to disk.
#' @param seed Randomization seed to use to choose additional sample signatures.
#'
#' @importFrom openxlsx read.xlsx
#' @import grDevices
#' @importFrom moments kurtosis
#'
#' @return No output.
#' @export
signatureDistributionPlot <- function(sigset="bhrr",
                              dataset="ph1_100normal_gene",
                              method="fc",
                              nullset = "ph1_100normal_gene_RAND125",
                              perc = .95,
                              fdr = .25,
                              comparetype = "Null",
                              samplepaths = c("HALLMARK_ESTROGEN_RESPONSE_EARLY", "DUTERTRE_ESTRADIOL_RESPONSE_6HR_UP",
                                              "HALLMARK_CHOLESTEROL_HOMEOSTASIS", "Vitamin A and carotenoid metabolism",
                                              "Cytochrome P450 signature", "HALLMARK_ANDROGEN_RESPONSE"),
                              to.file = T,
                              seed = 12345) {
  set.seed(seed)

  if(comparetype == "Null"){
    #load nullset scores
    file <- paste0("../output/signature_score_summary/signaturescoremat_",sigset,"_",nullset, "_",method,".RData")
    load(file)
    randoms = signaturescoremat

    #fill in sigset column. This might be deprecated.
    pway.annotated = read.xlsx("../input/all_signatures_annotated.xlsx")
    signaturescoremat$sigset = pway.annotated$sigset[match(signaturescoremat$signature, pway.annotated$signature)]
    randoms$sigset = pway.annotated$sigset[match(randoms$signature, pway.annotated$signature)]
  }

  #load dataset scores
  file <- paste0("../output/signature_score_summary/signaturescoremat_",sigset,"_",dataset,"_",method,".RData")
  load(file)

  #open pdf
  dir.create("../output/distplots/", showWarnings = FALSE)
  fname = paste0("../output/distplots/DISTPLOT_",sigset,"_",dataset,"_",method,"_",comparetype,".pdf")
  if(to.file) pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)

  #get random 10 extra signatures and add to user-specified signatures
  samplepaths = c(samplepaths,sample(unique(signaturescoremat$signature),10, replace = T))

  par(mfrow = c(2,1))
  cols = brewer.pal(6, "Set1")

  #highscores are the scores for actual data, comparescores are the null or low conc data
  highscores = signaturescoremat$signature_score
  if(comparetype == "Low Conc") comparescores = signaturescoremat$signature_score[signaturescoremat$conc < .2]
  if(comparetype == "Null") comparescores = randoms$signature_score

  #subtract comparescores median
  highscores = highscores - median(comparescores)
  comparescores = comparescores - median(comparescores)

  #compute pvals and fdr based on absolute scores
  myecdf = ecdf(abs(comparescores))
  pvals = 1-myecdf(abs(highscores))
  apvals = p.adjust(pvals, method = "BH")
  balfdrscore = min(abs(highscores)[apvals < fdr])

  #get densities for both
  maxscore = max(abs(highscores), abs(comparescores))
  dlen = 1024
  sdens = density(highscores, n = 1024, from = -maxscore, to = maxscore)
  ndens =  density(comparescores, n = 1024, from = -maxscore, to = maxscore)

  #plot densities
  ylims = c(0,max(c(sdens$y, ndens$y) ) )
  plot(sdens, xlab = "Pathway Score Probability Distribution",
       main = paste("Method:", method),ylim = ylims)
  points(ndens, col = cols[1], type = "l")

  #plot pval cutoff
  lines(rep(quantile(abs(comparescores),perc), 2),ylims, col = cols[2])
  lines(rep(-quantile(abs(comparescores),perc), 2),ylims, col = cols[2])

  #plot fdr cutoff
  if(is.finite(balfdrscore)){
    lines(rep(balfdrscore, 2),ylims, col = cols[5])
    lines(rep(-balfdrscore, 2),ylims, col = cols[5])
  }

  #get kurtosis for each and put in top left legend for densities
  label1 = paste0(comparetype, "(kurtosis = ", round(kurtosis(comparescores),1), ")")
  label2 = paste0("All Conc (kurtosis = ", round(kurtosis(highscores),1), ")")
  legend(x = "topleft", legend = c(label1, label2), col = c(cols[1],"black"), lwd = 2)

  #show top right legend for cutoffs
  if(is.finite(balfdrscore)){
    legend(x = "topright", legend = c(paste0("P-Value: ", 1-perc),
                                      paste0("FDR: ", fdr)), col = cols[c(2,5)], lty = 1)
  } else {
    legend(x = "topright", legend = c(paste0("P-Value: ", 1-perc),
                                      paste0("FDR: ", fdr, " (infinite)")), col = cols[c(2,5)], lty = 1)
  }

  #now repeat the above for all the samplepaths
  for(samplepath in samplepaths){
    highscores = signaturescoremat$signature_score[signaturescoremat$signature == samplepath]
    if(comparetype == "Low Conc") comparescores = signaturescoremat$signature_score[signaturescoremat$signature == samplepath & signaturescoremat$conc < .2]
    if(comparetype == "Null") comparescores = randoms$signature_score[randoms$signature == samplepath]

    highscores = highscores - median(comparescores)
    comparescores = comparescores - median(comparescores)

    myecdf = ecdf(abs(comparescores))
    pvals = 1-myecdf(abs(highscores))
    apvals = p.adjust(pvals, method = "BH")
    balfdrscore = min(abs(highscores)[apvals < fdr])

    maxscore = max(abs(highscores), abs(comparescores))
    dlen = 1024
    sdens = density(highscores, n = 1024, from = -maxscore, to = maxscore)
    ndens =  density(comparescores, n = 1024, from = -maxscore, to = maxscore)

    ylims = c(0,max(c(sdens$y, ndens$y) ) )
    plot(sdens, xlab = paste(samplepath,"Probability Distribution"),
         main = paste("Method:", method),ylim = ylims)
    points(ndens, col = cols[1], type = "l")

    label1 = paste0(comparetype)# (kurtosis = ", round(kurtosis(lowscores),1), ")")
    label2 = paste0("All Conc")# (kurtosis = ", round(kurtosis(highscores),1), ")")

    lines(rep(quantile(abs(comparescores),perc), 2),ylims, col = cols[2])
    lines(rep(-quantile(abs(comparescores),perc), 2),ylims, col = cols[2])

    if(is.finite(balfdrscore)){
      lines(rep(balfdrscore, 2),ylims, col = cols[5])
      lines(rep(-balfdrscore, 2),ylims, col = cols[5])
    }

    legend(x = "topleft", legend = c(label1, label2), col = c(cols[1],"black"), lwd = 2)

    if(is.finite(balfdrscore)){
      legend(x = "topright", legend = c(paste0("P-Value: ", 1-perc),
                                        paste0("FDR: ", fdr)), col = cols[c(2,5)], lty = 1)
    } else {
      legend(x = "topright", legend = c(paste0("P-Value: ", 1-perc),
                                        paste0("FDR: ", fdr, " (infinite)")), col = cols[c(2,5)], lty = 1)
    }
  }

  if(to.file) dev.off()
}
