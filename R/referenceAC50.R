#' Reference AC50 Plot
#'
#' Scatter plot and accuracy statistics of signatures vs. reference values.
#'
#' Saves a plot to disk. Plot is a scatter plot of actual values (based on
#' ER model) vs. predicted values (using some given signatures). For discrete
#' hitcalls, only true positive are plotted and colors indicate model used.
#' Continuous hitcalls plots all positives with colors indicating the hitcall.
#' Other statistics assume that all chemicals that are not positives (defined
#' by AUC >= .1) are negatives, so care must be taken not to include chemicals
#' with borderline activity in the dataset. RMSE is only shown for true
#' positives. Continuous hitcalls weights all statistics by the hitcall.
#' oldpval should be >= newpvals when using discrete hitcalls.
#'
#' @param method Pathway scoring method name.
#' @param dataset Data set name.
#' @param pathset Pathway set name.
#' @param nullset Null data set name.
#' @param newpvals Vector of p-values to make plots for.
#' @param oldpval P-value used when running signatureConcResp.
#' @param nametag Additional file identifier added during signatureConcResp.
#' @param conthits Set conthits = T when using continuous hits.
#' @param pathclass Some pre-defined sets of signatures to plot and run statistics
#'   on. "ER" is a group of ER signatures, "AR" is a group of AR signatures, and "DUT"
#'   is just the DUTERTRE_ESTRADIOL_RESPONSE_6HR_UP signature.
#' @param aucclass Which type of reference value to compare against. "erac50" uses
#'   the pseudo.AC50.median, "bmd" uses the pseudo.ACB.median, "AR" uses the
#'   maximum AR AUC, and "ER" uses the maximum ER AUC. AR, ER, and bmd might no
#'   longer function correctly.
#'
#' @import openxlsx
#' @import RColorBrewer
#'
#' @return No output.
#' @export
referenceAC50 = function(method = "fc", dataset = "user_wneg", pathset = "bhrr",
                         nullset = "user_wneg_RAND125",
                         newpvals = c(.2,.1,.05,.01,.005,.001), oldpval = .2,
                         nametag = NULL, conthits = F, pathclass = "DUT", aucclass = "erac50"){

  #preset signature classes
  if(pathclass == "ER") {
    keeppaths = c("DUTERTRE_ESTRADIOL_RESPONSE_6HR_UP", "HALLMARK_ESTROGEN_RESPONSE_EARLY",
                  "HALLMARK_ESTROGEN_RESPONSE_LATE", "RYAN_ESTROGEN_RECEPTOR_ALPHA_UP")
  } else if(pathclass == "AR") {
    keeppaths = c("Androgen receptor regulation of biosynthesis and transcription",
                  "Androgen receptor proteolysis and transcription regulation",
                  "HALLMARK_ANDROGEN_RESPONSE", "WANG_RESPONSE_TO_ANDROGEN_UP")
  } else if(pathclass == "DUT"){
    keeppaths = c("DUTERTRE_ESTRADIOL_RESPONSE_6HR_UP")
  } else {stop("Only pathclasses, DUT, ER and AR are defined.")}

  #read in er/ar scores for comparison
  erscores = read.xlsx("../extra/S2 ER SuperMatrix 2015-03-24.xlsx")
  arscores = read.xlsx("../extra/AR_supermatrix_Kleinstreuer_etal_2017.xlsx")

  #usescores is desired comparison metric
  keepcols = c("CODE", "CASRN", "Name", "AUC.Agonist", "AUC.Antagonist")
  if(aucclass == "ER") usescores= erscores[,keepcols]
  if(aucclass == "AR") usescores= arscores[,keepcols]
  if(aucclass == "erac50") usescores= erscores[,c(keepcols,"pseudo.AC50.median")]
  if(aucclass == "bmd") usescores= erscores[,c(keepcols,"pseudo.ACB.median")]

  #nametag formatting
  if(!is.null(nametag)) nametag = paste0("_", nametag)

  #load signature_CR
  file <- paste0("../output/signature_conc_resp_summary/PATHWAY_CR_",pathset,"_",dataset,"_",method ,"_",oldpval,nametag,
                 ".RData")
  load(file)
  PATHWAY_CR = PATHWAY_CR[PATHWAY_CR$signature %in% keeppaths,]

  #open output pdf
  dir.create("../output/refchems/", showWarnings = F)
  file = paste0("../output/refchems/refAC50_",pathset,"_",dataset,"_",method,"_",pathclass,"_",
              aucclass,nametag, ".pdf")
  pdf(file= file,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  par(mar = c(5.1,4.1,4.1,1.1),mfrow = c(3,2))

  #color by fitmethod
  methkey = c("poly1", "poly2", "pow", "exp2", "exp3", "exp4", "exp5", "hill", "gnls")
  methcols = c( brewer.pal(5, "Blues"), brewer.pal(3,"Reds"), "green")
  names(methcols) = methkey

  #get max ag/antag for ers/ar AUCs >= cutoff of .1
  erscores$maxauc = pmax(erscores$AUC.Agonist, erscores$AUC.Antagonist)
  erantags = erscores$CASRN[erscores$AUC.Antagonist > erscores$AUC.Agonist]
  ers = erscores$CASRN[erscores$maxauc >= .1]
  arscores$maxauc = pmax(arscores$AUC.Agonist, arscores$AUC.Antagonist)
  arantags = arscores$CASRN[arscores$AUC.Antagonist > arscores$AUC.Agonist]
  ars = arscores$CASRN[arscores$maxauc >= .1]

  #cycle through new pvalues, create one plot per
  for(newpval in newpvals){
    newpcr = PATHWAY_CR

    #recompute hitcalls for new pvalue
    pvalkey = getpvalcutoff(pathset, nullset = nullset, method = method, pvals = newpval)
    newcutoff = pvalkey$cutoff[match(newpcr$signature, pvalkey$signature)]
    if(conthits) newpcr$hitcall = hitcont(newpcr, newcutoff = newcutoff) else {
      newpcr$hitcall = hitlogic(newpcr, newcutoff = newcutoff)
    }

    #use only chemicals present in both data and reference
    keepchems = intersect(usescores$CASRN, newpcr$casrn)
    usescores = usescores[usescores$CASRN %in% keepchems,]
    newpcr = newpcr[newpcr$casrn %in% keepchems,]

    #get max auc (for when aucclass == ER or AR), "maxer" is overwritten if not
    newpcr$erags = usescores$AUC.Agonist[match(newpcr$casrn ,usescores$CASRN)]
    newpcr$erants = usescores$AUC.Antagonist[match(newpcr$casrn ,usescores$CASRN)]
    newpcr$maxer = pmax(newpcr$erags, newpcr$erants)

    #fill in reference values for other aucclasses
    if(aucclass == "erac50"){
      newpcr$erac50 = usescores$pseudo.AC50.median[match(newpcr$casrn ,usescores$CASRN)]
      newpcr$maxer = -log10(newpcr$erac50)
      yname = "-log10(Pseudo AC50 Median)"
      xname = "Predicted -log10(AC50)"
    } else if(aucclass == "bmd"){
      newpcr$erac50 = usescores$pseudo.ACB.median[match(newpcr$casrn ,usescores$CASRN)]
      newpcr$maxer = -log10(newpcr$erac50)
      yname = "-log10(Pseudo BMD Median)"
      xname = "Predicted -log10(BMD)"
    } else {
      yname = paste0("Max AUC ",aucclass," Agonist/Antagonist")
      xname = "Predicted -log10(AC50)"
    }

    if(conthits == T) newhits = -1 else newhits = .5

    #the newhit value decided the hitcall threshhold for points to plot; it's largely deprecated
    for(newhit in newhits){
      #fill in predicted values to compare
      par(xpd = F)
      if(aucclass == "bmd") newpcr$nlac50 = -log10(newpcr$bmd) else newpcr$nlac50 = -log10(newpcr$ac50)
      newpcr$nlac50[newpcr$hitcall < newhit] = NA

      pathsymbs = 21:24

      #set plot ranges
      yrange = range(newpcr$maxer, na.rm = T)
      yrange[1] = min(0,yrange[1])
      yrange[2] = max(1,yrange[2])
      xrange = range(newpcr$nlac50, na.rm = T)
      xrange[1] = min(-2,xrange[1])
      xrange[2] = max(2.5,xrange[2])
      if(aucclass == "erac50" || aucclass == "bmd") xrange = yrange = c(-2,3)

      #P/N by signature for discrete hitcalls
      P = sum(newpcr$casrn %in% ers)
      TP =  sum(newpcr$casrn %in% ers & (!is.na(newpcr$nlac50)))
      N = sum(!newpcr$casrn %in% ers)
      TN = sum((is.na(newpcr$nlac50)) & (!newpcr$casrn %in% ers))
      RMSE = sqrt(mean( ((newpcr$nlac50 - newpcr$maxer)[newpcr$casrn %in% ers])^2 , na.rm = T))
      Rsq = R2(newpcr$maxer[newpcr$casrn %in% ers], newpcr$nlac50[newpcr$casrn %in% ers])

      #P/N for continuous hitcalls
      if(conthits){
        P = sum(newpcr$casrn %in% ers)
        TP =  sum(as.numeric(newpcr$casrn %in% ers)*newpcr$hitcall)
        FP =  sum(as.numeric(!newpcr$casrn %in% ers)*newpcr$hitcall)
        N = sum(!newpcr$casrn %in% ers)
        TN =  sum(as.numeric(!newpcr$casrn %in% ers)*(1-newpcr$hitcall))
        FN =  sum(as.numeric(newpcr$casrn %in% ers)*(1-newpcr$hitcall))

        #weighted RMSE and R2 by hitcall
        pospcr = newpcr[newpcr$casrn %in% ers & newpcr$hitcall > 0,]
        weight = pospcr$hitcall
        RMSE = sqrt(sum(((pospcr$nlac50 - pospcr$maxer)^2)*weight)/sum(weight) )

        totvar = sum((pospcr$maxer-mean(pospcr$maxer))^2*weight)/sum(weight)
        Rsq = (1 - RMSE^2/totvar)

        #continuous hitcall plot setup
        newpcr$contcol = "black"
        allcols = brewer.pal(9,"Spectral")
        hitbnds = c(1,.99,.9,.5,.1,.01,.001,.0001,.00001)
        newpcr$contcol = sapply(newpcr$hitcall, function(x){allcols[max(which(x < hitbnds))]})
        plot(newpcr$nlac50, newpcr$maxer, col = "black",
             bg = newpcr$contcol, pch = pathsymbs[as.factor(newpcr$signature)],
             xlab = xname, ylab = yname, ylim = yrange, xlim = xrange)
      } else {
        #discrete hitcall plot setup
        plot(newpcr$nlac50, newpcr$maxer, col = "black",
           bg = methcols[newpcr$fit_method], pch = pathsymbs[as.factor(newpcr$signature)],
           xlab = xname, ylab = yname, ylim = yrange, xlim = xrange)
      }

      #title, including statistics
      headtitle = "P-Value: "
      if(aucclass == "erac50" | aucclass == "bmd"){
        title(paste0(headtitle, newpval, ", BA: ", round(mean(c(TP/P, TN/N)),2),
                     "\nTPR: ",round(TP/P,2), ", TNR: ", round(TN/N,2),
                     "\nR2: ", round(Rsq,2),", RMSE: ", round(RMSE,2) ),
              adj = 0, line = .5, cex.sub = .7)
      } else {
        title(paste0(headtitle, newpval, ", Hit:", newhit,"\nTPR: ",round(TP/P,2), " ,TNR: ", round(TN/N,2), "\nBA: ",
              round(mean(c(TP/P, TN/N)),2) ), adj = 0, line = .5, cex.sub = .7)
      }

      #y=x reference line with 1 log unit dotted bounds
      lines(xrange,yrange)
      if(aucclass == "erac50" || aucclass == "bmd") {
        onesd = 1
        lines(xrange, yrange-onesd, lty = 2)
        lines(xrange, yrange + onesd, lty = 2)
      }

      #legends
      par(xpd = T)
      legend("topright", inset = c(0,-.23), legend = tolower(strtrim(keeppaths, 40)), col = "black",
           pch = pathsymbs, cex = .7)
      if(conthits) {
        colkey = paste0("<" , hitbnds[1:9])
        legend("bottomright", inset = c(0,0), legend = colkey, fill = allcols, cex = .7)
      }else {
        legend("bottomright", inset = c(0,0), legend = methkey, fill = methcols, cex = .7)
      }

    }
  }
  dev.off()


}
