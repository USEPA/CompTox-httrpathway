#--------------------------------------------------------------------------------------
#' Calculate several statistics on a 2 x 2 matrix
#'
#' @param tp number of true positives
#' @param fp number of false positives
#' @param fn number of false negatives
#' @param tn number of true negatives
#' @param do.p if TRUE, calcualte an exact p-value
#' @param rowname if not NA, adda column to the output with this rowname
#'
#' Returns:
#'   a list of the results
#'   a: TP
#'   b: FP
#'   c: FN
#'   d: TN
#'   sens: sensitivity
#'   spec: specificity
#'   ba: Balanced Accuracy
#'   accuracy: Accuracy
#'   relative.risk: Relative Risk
#'   odds.ratio: Odds Ratio
#'   or.ci.lwr: lower confidence interval of the Odds Ratio
#'   or.ci.upr: upper confidence interval of the Odds Ratio
#'   ppv: Positive Predictive Value
#'   npv: Negative Predictive Value
#'   p.value: Chi-squared p-value
#'   F1: 2TP/(2TP+FP+FN)
#'
#'   sval: All of the results as a tab-delimited string
#'   title: the title of the results as a tab-delimited string
#'  mat: The results as a 1-row data frame
#'  @export
#--------------------------------------------------------------------------------------
TxT <<- function(tp,fp,fn,tn,do.p=TRUE,rowname=NA) {
  sens <- tp/(tp+fn)
  spec <- tn/(tn+fp)
  ppv  <- tp/(tp+fp)
  npv  <- tn/(tn+fn)
  F1 <- 2*tp/(2*tp+fp+fn)

  relative.risk <- (tp/(tp+fp)) / (fn/(tn+fn))
  odds.ratio <- (tp*tn)/(fp*fn)

  accuracy <- (tp+tn)/(tp+tn+fp+fn)
  x<-matrix(data=NA, nrow=2, ncol=2)
  x[1,1]<-tp
  x[1,2]<-fp
  x[2,1]<-fn
  x[2,2]<-tn
  if(is.infinite(relative.risk))  relative.risk <- 1000000
  if(is.infinite(odds.ratio))  odds.ratio <- 1000000
  if(is.infinite(sens))  sens <- 0
  if(is.infinite(spec))  spec <- 0
  if(is.nan(relative.risk))  relative.risk <- 1
  if(is.nan(odds.ratio))  odds.ratio <- 1
  if(is.nan(sens))  sens <- 0
  if(is.nan(spec))  spec <- 0
  if(is.na(sens))  sens <- 0
  if(is.na(spec))  spec <- 0
  or.ci.lwr <- 0
  or.ci.upr <- 1000000
  if(odds.ratio<1000000 && tp>0 && tn>0 && fp>0 && fn>0) {
    lnor <- log(odds.ratio)
    selnor <- sqrt(1/tp+1/tn+1/fp+1/fn)
    ln.upr <- lnor+1.96*selnor
    ln.lwr <- lnor-1.96*selnor
    or.ci.lwr <- exp(ln.lwr)
    or.ci.upr <- exp(ln.upr)
  }
  p.value<-1
  if(do.p==TRUE) {
    if(sens>0 && spec>0) {
      c<-fisher.test(x)
      p.value <- c$p.value
    }
  }
  ba <- 0.5*(sens+spec)

  sval<-paste(tp,"\t",fp,"\t",fn,"\t",tn,"\t",
              format(sens,digits=3),"\t",
              format(spec,digits=3),"\t",
              format(ba,digits=3),"\t",
              format(accuracy,digits=3),"\t",
              format(relative.risk,digits=3),"\t",
              format(odds.ratio,digits=3),"\t",
              format(or.ci.lwr,digits=3),"\t",
              format(or.ci.upr,digits=3),"\t",
              format(ppv,digits=3),"\t",
              format(npv,digits=3),"\t",
              format(F1,digits=3),"\t",
              format(p.value,digits=3),sep="")
  title<-paste("TP\tFP\tFN\tTN\tSens\tSpec\tBA\tAcrcy\tRelRsk\tOR\tCI.OR.LWR\tCI.OR.UPR\tPPV\tNPV\tF1\tp.value",sep="")

  name.list <- c("TP","FP","FN","TN","Sens","Spec","BA","Acrcy","RelRsk","OR","CI.OR.LWR","CI.OR.UPR","PPV","NPV","F1","p.value")
  if(!is.na(rowname)) name.list <- c("rowname",name.list)
  mat <- as.data.frame(matrix(nrow=1,ncol=length(name.list)))
  names(mat) <- name.list
  mat[1,"TP"] <- tp
  mat[1,"FP"] <- fp
  mat[1,"FN"] <- fn
  mat[1,"TN"] <- tn
  mat[1,"Sens"] <- sens
  mat[1,"Spec"] <- spec
  mat[1,"BA"] <- ba
  mat[1,"Acrcy"] <- accuracy
  mat[1,"RelRsk"] <- relative.risk
  mat[1,"OR"] <- odds.ratio
  mat[1,"CI.OR.LWR"] <- or.ci.lwr
  mat[1,"CI.OR.UPR"] <- or.ci.upr
  mat[1,"PPV"] <- ppv
  mat[1,"NPV"] <- npv
  mat[1,"F1"] <- F1
  mat[1,"p.value"] <- p.value
  if(!is.na(rowname)) mat[1,"rowname"] <- rowname

  r <- list(a=tp,b=fp,c=fn,d=tn,sens=sens,spec=spec,ba=ba,accuracy=accuracy,relative.risk=relative.risk,odds.ratio=odds.ratio,or.ci.lwr=or.ci.lwr,or.ci.upr=or.ci.upr,ppv=ppv,npv=npv,F1=F1,p.value=p.value,sval=sval,title=title,mat=mat)
  r
}
