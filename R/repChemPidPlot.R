#' Replicate Chemical PID Plot
#'
#' Generates plots and statistics for replicate chemicals' probe IDs.
#'
#' This function is desigend to work with runAllRepChemPidCR, so the dataset names
#' are hard-coded. This function may take some time to run. Concentration response
#' should have been run using continuous hitcalls.
#'
#' @param oldpval P-value used to generate GENE_CR's.
#' @param nametag Optional descriptor in filename.
#' @param mc.cores Number of cores to use.
#'
#'@import openxlsx
#'@import RColorBrewer
#'@import parallel
#'
#' @return No output.
#' @export
repChemPidPlot = function(oldpval = .05, nametag = "conthits", mc.cores = 3){

  #NOTE: see repChemPathwayPlot() for more comments; code is nearly identical.
  if(!is.null(nametag)) nametag = paste0("_", nametag)

  floors = c("5","10")
  pes = c("0","1")
  methods = c("none", "normal", "normalold", "apeglm", "ashr")
  combos = expand.grid(floors,pes,methods, stringsAsFactors = F)
  datanames = apply(combos,1,function(x){paste0(c(x,"_pid"),collapse = "")})

  erscores = read.xlsx("S2 ER SuperMatrix 2015-03-24.xlsx")
  pos = erscores$CASRN[erscores$AUC.Antagonist >= .1 | erscores$AUC.Agonist >=.1]
  neg = erscores$CASRN[erscores$pseudo.AC50.median == 1000000]
  load("../input/processed_signature_data/PATHWAY_LIST_bhrr.RData")
  genes = signature_data[["DUTERTRE_ESTRADIOL_RESPONSE_6HR_UP"]]

  cl = makePSOCKcluster(mc.cores)
  clusterExport(cl, c("RMSE", "WRMSE", "auc"))

  ###########################################Data Load and Process################################################
  finalout = parLapply(cl = cl, X=datanames, fun= function(dataset, oldpval, nametag, pos, neg, genes){

    file <- paste0("../output/gene_conc_resp_summary/GENE_CR_ph1_",dataset,"_",oldpval,nametag,
                   ".RData")
    load(file)

    ph1 = GENE_CR
    #ph1 has a duplicate PFOS; easiest just to ignore one of them
    ph1 = ph1[ph1$sample_id != "TP0001718N24",]

    file <- paste0("../output/gene_conc_resp_summary/GENE_CR_pilot_",dataset,"_",oldpval,nametag,
                   ".RData")
    load(file)

    pilot = GENE_CR

    pilotid = paste0(pilot$dtxsid, "_", pilot$gene)
    ph1id = paste0(ph1$dtxsid, "_", ph1$gene)
    commonid = intersect(pilotid, ph1id)
    pilot = pilot[pilotid %in% commonid,]
    ph1 = ph1[ph1id %in% commonid,]

    #order, so both are one-to-one
    pilot = pilot[order(pilot$dtxsid, pilot$gene),]
    ph1 = ph1[order(ph1$dtxsid, ph1$gene),]

    hitcor = cor(ph1$hitcall, pilot$hitcall)
    hitcor2 = cor(ph1$hitcall, pilot$hitcall, method = "spearman")

    quants = quantile(c(pilot$hitcall,pilot$hitcall), 1:100/100)
    quants = unique(quants)

    hitmatch = sapply(quants, function(x){ sum((ph1$hitcall >= x) == (pilot$hitcall >= x)) })
    hitnums = sapply(quants, function(x){ sum(ph1$hitcall >= x) + sum(pilot$hitcall >= x) })
    hitfracs = hitnums/(nrow(pilot) + nrow(ph1))
    matchfracs = hitmatch/nrow(pilot)

    hitboth = sapply(quants, function(x){ sum((ph1$hitcall >= x) & (pilot$hitcall >= x)) })
    bothfrac = hitboth/nrow(pilot)

    qrmseac50 = sapply(quants, function(x){
      return(WRMSE(rep(log10(pilot$ac50),2), rep(log10(ph1$ac50),2), w = (c(pilot$hitcall,ph1$hitcall) > x)))
    })
    names(qrmseac50) = quants

    qrmsebmd = sapply(quants, function(x){
      return(WRMSE(rep(log10(pilot$bmd),2), rep(log10(ph1$bmd),2), w = (c(pilot$hitcall,ph1$hitcall) > x)))
    })

    names(qrmsebmd) = quants

    qrmseac50both = sapply(quants, function(x){
      return(WRMSE(log10(pilot$ac50), log10(ph1$ac50), w = ((pilot$hitcall > x) & (ph1$hitcall > x)) ))
    })
    names(qrmseac50both) = quants

    qrmsebmdboth = sapply(quants, function(x){
      return(WRMSE(log10(pilot$bmd), log10(ph1$bmd), w = ((pilot$hitcall > x) & (ph1$hitcall > x)) ))
    })
    names(qrmsebmdboth) = quants

    #Phase 1 Estrogen Detection
    erph1 = ph1[ph1$casrn %in% c(pos,neg),]
    erph1$gene = sub("_.*","",erph1$gene)
    erph1 = erph1[erph1$gene %in% genes,]
    P =  sum(erph1$casrn %in% pos)
    N =  sum(erph1$casrn %in% neg)
    TPS = sapply(quants, function(i){ sum(erph1$hitcall[erph1$casrn %in% pos] >= i) })
    FPS = sapply(quants, function(i){ sum(erph1$hitcall[erph1$casrn %in% neg] >= i) })
    TNS = sapply(quants, function(i){ sum(erph1$hitcall[erph1$casrn %in% neg] < i)  })
    FNS = sapply(quants, function(i){ sum(erph1$hitcall[erph1$casrn %in% pos] < i)  })
    BAS = (TPS/P + TNS/N)/2
    tpr = TPS/P
    fpr = FPS/N
    AUC = auc(tpr, fpr)
    ph1erstats = list(tpr = tpr, fpr = fpr, ba = BAS, AUC = AUC)

    #Pilot Estrogen Detection
    erpilot = pilot[pilot$casrn %in% c(pos,neg),]
    erpilot$gene = sub("_.*","",erpilot$gene)
    erpilot = erpilot[erpilot$gene %in% genes,]
    P =  sum(erpilot$casrn %in% pos)
    N =  sum(erpilot$casrn %in% neg)
    TPS = sapply(quants, function(i){ sum(erpilot$hitcall[erpilot$casrn %in% pos] >= i) })
    FPS = sapply(quants, function(i){ sum(erpilot$hitcall[erpilot$casrn %in% neg] >= i) })
    TNS = sapply(quants, function(i){ sum(erpilot$hitcall[erpilot$casrn %in% neg] < i)  })
    FNS = sapply(quants, function(i){ sum(erpilot$hitcall[erpilot$casrn %in% pos] < i)  })
    BAS = (TPS/P + TNS/N)/2
    tpr = TPS/P
    fpr = FPS/N
    AUC = auc(tpr, fpr)
    piloterstats = list(tpr = tpr, fpr = fpr, ba = BAS, AUC = AUC)

    return(list( hitcor = hitcor, hitcor2 = hitcor2, hitfrac = hitfracs, matchfrac = matchfracs, bothfrac= bothfrac, quants = quants,
                 qrmseac50 = qrmseac50, qrmsebmd = qrmsebmd, qrmseac50both = qrmseac50both, qrmsebmdboth = qrmsebmdboth,
                 ph1erstats = ph1erstats, piloterstats = piloterstats))


  }, oldpval = oldpval, nametag = nametag, pos = pos, neg = neg, genes = genes)

  stopCluster(cl)
  names(finalout) = datanames

  ###############################################Data Extraction####################################################

  quants = sapply(finalout, function(totalout){ return(totalout$quants)})
  hitfrac = sapply(finalout, function(totalout){ return(totalout$hitfrac)})
  hitcor = sapply(finalout, function(totalout){ return(totalout$hitcor)})
  hitcor2 = sapply(finalout, function(totalout){ return(totalout$hitcor2)})
  matchfrac = sapply(finalout, function(totalout){ return(totalout$matchfrac)})
  qrmseac50 = sapply(finalout, function(totalout){ return(totalout$qrmseac50)})
  qrmsebmd = sapply(finalout, function(totalout){ return(totalout$qrmsebmd)})

  bothfrac = sapply(finalout, function(totalout){ return(totalout$bothfrac)})
  qrmseac50both = sapply(finalout, function(totalout){ return(totalout$qrmseac50both)})
  qrmsebmdboth = sapply(finalout, function(totalout){ return(totalout$qrmsebmdboth)})

  aucph1 = sapply(finalout, function(totalout){ return(totalout$ph1erstats$AUC)})
  aucpilot = sapply(finalout, function(totalout){ return(totalout$piloterstats$AUC)})
  baph1 = sapply(finalout, function(totalout){ return(totalout$ph1erstats$ba)})
  maxbaph1 = sapply(baph1,max, na.rm = T)
  bapilot = sapply(finalout, function(totalout){ return(totalout$piloterstats$ba)})
  maxbapilot = sapply(bapilot,max, na.rm = T)

  tprph1 = sapply(finalout, function(totalout){ return(totalout$ph1erstats$tpr)})
  fprph1 = sapply(finalout, function(totalout){ return(totalout$ph1erstats$fpr)})

  tprpilot = sapply(finalout, function(totalout){ return(totalout$piloterstats$tpr)})
  fprpilot = sapply(finalout, function(totalout){ return(totalout$piloterstats$fpr)})

  #################################################Plotting#############################################################

  file = "../output/repchem/pidplots.pdf"
  pdf(file= file,width=8.5,height=11,pointsize=12,bg="white",paper="letter",pagecentre=T)
  par(mar = c(4.1,4.1,2.1,1.1),mfrow = c(2,1))

  cols = brewer.pal(10,"Paired")
  cols = rep(cols, each = 2)
  types = rep(1:2,10)

  labs = sub("_pid", "", datanames)
  labs = sub("50", "5_0_", labs)
  labs = sub("100", "10_0_", labs)
  labs = sub("51", "5_1_", labs)
  labs = sub("101", "10_1_", labs)

  #hitfrac vs. hitcall
  plot(0,type = "n", ylim = c(0,1), xlim = c(0,1.3), ylab = "Combined Hit Fraction", xlab = "Hitcall")
  silentout = mapply(function(x,y,col,type){points(x,y,type = "l", lty = type, col = col, lwd = 2)},
                     y = hitfrac, x = quants, col = cols, type = types)
  legend("right", legend = labs, lty = types, col = cols, cex = .85, lwd = 2)

  #matchfrac vs. hitfrac
  plot(0,type = "n", ylim = c(.5,1), xlim = c(0,1.3), ylab = "Fraction of Hitcalls Matching", xlab = "Combined Hit Fraction")
  silentout = mapply(function(x,y,col,type){points(x,y,type = "l", lty = type, col = col, lwd = 2)},
                     y = matchfrac, x = hitfrac, col = cols, type = types)
  x = 0:100/100
  y = x^2 + (1-x)^2   #Perfectly random distribution
  points(x,y,type = "l", lwd = 2)
  legend("right", legend = labs, lty = types, col = cols, cex = .85, lwd = 2)

  #bothfrac vs. hitfrac
  plot(0,type = "n", ylim = c(0,1), xlim = c(0,1.3), ylab = "Fraction of Pathways Hitting Both Replicates",
       xlab = "Combined Hit Fraction")
  silentout = mapply(function(x,y,col,type){points(x,y,type = "l", lty = type, col = col, lwd = 2)},
                     y = bothfrac, x = hitfrac, col = cols, type = types)
  x = 0:100/100
  y = x^2   #Perfectly random distribution
  points(x,y,type = "l", lwd = 2)
  legend("right", legend = labs, lty = types, col = cols, cex = .85, lwd = 2)

  #AC50 replication
  plot(0,type = "n", ylim = c(0,2.5), xlim = c(0,1.3), ylab = "AC50 Replication RMSE", xlab = "Combined Hit Fraction")
  silentout = mapply(function(x,y,col,type){points(x,y,type = "l", lty = type, col = col, lwd = 2)},
         y = qrmseac50, x = hitfrac, col = cols, type = types)
  legend("right", legend = labs, lty = types, col = cols, cex = .85, lwd = 2)

  #BMD replication
  plot(0,type = "n", ylim = c(0.5,2.0), xlim = c(0,1.3), ylab = "BMD Replication RMSE", xlab = "Combined Hit Fraction")
  silentout = mapply(function(x,y,col,type){points(x,y,type = "l", lty = type, col = col, lwd = 2)},
                     y = qrmsebmd, x = hitfrac, col = cols, type = types)
  legend("right", legend = labs, lty = types, col = cols, cex = .85, lwd = 2)

  #AC50 both replication
  plot(0,type = "n", ylim = c(0,2.5), xlim = c(0,1.3), ylab = "AC50 Replication RMSE (Both Hit)", xlab = "Combined Hit Fraction")
  silentout = mapply(function(x,y,col,type){points(x,y,type = "l", lty = type, col = col, lwd = 2)},
                     y = qrmseac50both, x = hitfrac, col = cols, type = types)
  legend("right", legend = labs, lty = types, col = cols, cex = .85, lwd = 2)

  #BMD both replication
  plot(0,type = "n", ylim = c(0.5,2.5), xlim = c(0,1.3), ylab = "BMD Replication RMSE (Both Hit)", xlab = "Combined Hit Fraction")
  silentout = mapply(function(x,y,col,type){points(x,y,type = "l", lty = type, col = col, lwd = 2)},
                     y = qrmsebmdboth, x = hitfrac, col = cols, type = types)
  legend("right", legend = labs, lty = types, col = cols, cex = .85, lwd = 2)

  # #ROC Combined
  # plot(0,type = "n", ylim = c(0,1), xlim = c(0,1.3), ylab = "True Positive Rate", xlab = "False Positive Rate",
  #      main = "ROC Curve - Pilot and Phase 1 Combined")
  # silentout = mapply(function(x,y,col,type){points(x,y,type = "l", lty = type, col = col, lwd = 2)},
  #                    y = tpr, x = fpr, col = cols, type = types)
  # legend("right", legend = labs, lty = types, col = cols, cex = .85, lwd = 2)
  # abline(0,1, lwd = 2)

  #ROC Phase 1
  plot(0,type = "n", ylim = c(0,1), xlim = c(0,1.3), ylab = "True Positive Rate", xlab = "False Positive Rate",
       main = "ROC Curve - Phase 1")
  silentout = mapply(function(x,y,col,type){points(x,y,type = "l", lty = type, col = col, lwd = 2)},
                     y = tprph1, x = fprph1, col = cols, type = types)
  abline(0,1, lwd = 2)
  legend("right", legend = labs, lty = types, col = cols, cex = .85, lwd = 2)


  #ROC Pilot
  plot(0,type = "n", ylim = c(0,1), xlim = c(0,1.3), ylab = "True Positive Rate", xlab = "False Positive Rate",
       main = "ROC Curve - Pilot")
  silentout = mapply(function(x,y,col,type){points(x,y,type = "l", lty = type, col = col, lwd = 2)},
                     y = tprpilot, x = fprpilot, col = cols, type = types)
  abline(0,1, lwd = 2)
  legend("right", legend = labs, lty = types, col = cols, cex = .85, lwd = 2)


  #Phase 1 BA vs.hitfrac
  plot(0,type = "n", ylim = c(.2,.6), xlim = c(0,1.3), ylab = "Phase 1 Balanced Accuracy", xlab = "Combined Hit Fraction")
  silentout = mapply(function(x,y,col,type){points(x,y,type = "l", lty = type, col = col, lwd = 2)},
                     y = baph1, x = hitfrac, col = cols, type = types)
  abline(.5,0, lwd = 2)
  legend("right", legend = labs, lty = types, col = cols, cex = .85, lwd = 2)

  #Pilot BA vs. hitfrac
  plot(0,type = "n", ylim = c(.4,.8), xlim = c(0,1.3), ylab = "Pilot Balanced Accuracy", xlab = "Combined Hit Fraction")
  silentout = mapply(function(x,y,col,type){points(x,y,type = "l", lty = type, col = col, lwd = 2)},
                     y = bapilot, x = hitfrac, col = cols, type = types)
  abline(.5,0, lwd = 2)
  legend("right", legend = labs, lty = types, col = cols, cex = .85, lwd = 2)

  #BAR PLOTS
  #hitcall correlation
  pidbar(hitcor, ylab = "Pearson Correlation of Hitcalls")
  pidbar(hitcor2, ylab = "Spearman Correlation of Hitcalls")

  #ER AUCs
  pidbar(aucph1, ylab = "AUC Phase 1 ERs")
  pidbar(aucpilot, ylab = "AUC Pilot ERs")
  pidbar(rbind(aucph1, aucpilot), ylab = "AUC Combined (Phase 1 Bottom)")

  #ER BAs
  pidbar(maxbaph1, ylab = "Phase 1 Maximum BA ERs")
  pidbar(maxbapilot, ylab = "Pilot Maximum BA ERs")
  pidbar(rbind(maxbaph1, maxbapilot), ylab = "Combined Maximum BA ERs (Phase 1 Bottom)")

  dev.off()

}

#' PID Bar Plot
#'
#' Specially formatted bar plot.
#'
#' This function is a helper for repChemPidPlot. It fiddles with the margins
#' and renames the labels so that they fit on the plot.
#'
#' @param x Named matrix or vector to pass to barplot.
#' @param ... Other options to pass to barplot.
#'
#' @return No output.
pidbar = function(x, ...){
  par(mar = c(7.6,4.1,1.1,1.1))
  if(is.matrix(x)) labs = colnames(x) else labs = names(x)
  labs = sub("_pid", "", labs)
  labs = sub("50", "5_0_", labs)
  labs = sub("100", "10_0_", labs)
  labs = sub("51", "5_1_", labs)
  labs = sub("101", "10_1_", labs)
  barplot(x, names.arg = labs, las = 2, cex.names = 1,...)
  par(mar = c(4.1,4.1,2.1,1.1))
}
