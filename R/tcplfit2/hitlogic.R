#' Hit Logic (Discrete)
#' 
#' Wrapper that computes discrete hitcalls for a provided PATHWAY_CR dataframe.
#'
#' @param indf Dataframe similar to PATHWAY_CR. Must contain "conc" and "resp"
#'   columns if xs and ys are not provided. Must contain "cutoff" and "bmad_factor"
#'   columns if newbmad is not NULL. Must contain "top" and "ac50" columns. "conc"
#'   and "resp" entries should be a single string with values separated by |.
#' @param newbmad (Deprecated) New number of bmads to use for the cutoff.
#' @param xs List of concentration vectors that can be provided for speed.
#' @param ys List of response vectors that can be provided for speed.
#' @param newcutoff Vector of new cutoff values to use. Length should be equal
#'   to rows in indf. 
#'
#' @return Vector of hitcalls with length equal to number of rows in indf.
#' @export
#'
#' @examples
#' conc = rep(".03|.1|.3|1|3|10|30|100",2)
#' resp = rep("0|0|.1|.1|.5|.5|1|1",2)
#' indf = data.frame(top = c(1,1), ac50 = c(3,4), conc = conc, resp = resp, 
#'   stringsAsFactors = FALSE)
#' hitlogic(indf, newcutoff = c(.8,1.2))
hitlogic = function(indf, newbmad = NULL, xs = NULL, ys = NULL, newcutoff = NULL){
  
  #extract cutoff from newbmad or newcutoff
  if(!is.null(newbmad)) cutoff = indf$cutoff/indf$bmad_factor*newbmad
  if(!is.null(newcutoff)) cutoff = newcutoff
  
  #reformat concs and resps, if necessary
  if(is.null(xs)){
    xs = strsplit(indf$conc, "\\|")
    xs = lapply(xs, as.numeric)
  }
  if(is.null(ys)){
    ys = strsplit(indf$resp, "\\|")
    ys = lapply(ys, as.numeric)
  }
  
  #run hitoginner for each row of indf
  hitcall = mapply(hitloginner, conc= xs, resp = ys, top = indf$top, cutoff = cutoff, ac50 = indf$ac50)
  
  return(hitcall)
}


#' Hit Logic Inner (Discrete)
#' 
#' Contains hit logic, called directly during CR fitting or later through "hitlogic".
#' 
#' The purpose of this function is to keep the actual hit rules in one
#' location so it can be called during CR fitting, and then again after the fact
#' for a variety of cutoffs. Curves fit with constant winning should have
#' top = NA, generating a miss.
#'
#' @param conc Vector of concentrations (No longer necessary).
#' @param resp Vector of responses.
#' @param top Model top.
#' @param cutoff Desired cutoff.
#' @param ac50 Model AC50 (No longer necessary).
#'
#' @return Outputs 1 for hit, 0 for miss.
#' @export
#'
#' @examples
#' hitloginner(resp = 1:8, top = 7, cutoff = 5) #hit
#' hitloginner(resp = 1:8, top = 7, cutoff = 7.5) #miss: top too low
#' hitloginner(resp = 1:8, top = 9, cutoff = 8.5) #miss: no response> cutoff
#' hitloginner(resp = 1:8, top = NA, cutoff = 5) #miss: no top (constant)
hitloginner = function(conc = NULL, resp, top, cutoff, ac50 = NULL){
  
  
  n_gt_cutoff = sum(abs(resp)>cutoff)
  
  #hitlogic - hit must have: at least one point above abs cutoff, 
  # a defined top (implying there is a winning non-constant model),
  #and an abs. top greater than the cutoff
  hitcall = 0
  if(n_gt_cutoff>0 && !is.na(top) && abs(top)>cutoff) hitcall <- 1
  
  return(hitcall)
}

#' Continuous Hitcalls
#' 
#' Wrapper that computes continuous hitcalls for a provided PATHWAY_CR dataframe.
#' 
#' indf parameter columns should be NA when not required by fit method. "conc" 
#' and "resp" entries should be a single string with values separated by |. 
#' Details on indf columns can be found in pathwayConcRespCore_pval.
#'
#' @param indf Dataframe similar to PATHWAY_CR. Must contain "conc" and "resp"
#'   columns if xs and ys are not provided. Must contain "top", "ac50", "er",
#'   "fit_method", "caikwt", and "mll" columns as well as columns for each
#'   model parameter. 
#' @param xs List of concentration vectors that can be provided for speed.
#' @param ys List of response vectors that can be provided for speed.
#' @param newcutoff Vector of new cutoff values to use. Length should be equal
#'   to rows in indf. 
#' @param mc.cores Number of cores to use for large dataframes.
#' 
#' @import future.apply
#' @import future
#'
#' @return Vector of hitcalls between 0 and 1 with length equal to indf row
#'   number.
#' @export
hitcont = function(indf, xs = NULL, ys = NULL, newcutoff, mc.cores = 1){
  
  # reformat concs and resps
  if(is.null(xs)){
    xs = strsplit(indf$conc, "\\|")
    xs = lapply(xs, as.numeric)
  }
  if(is.null(ys)){
    ys = strsplit(indf$resp, "\\|")
    ys = lapply(ys, as.numeric)
  }
  
  #correct parameter ordering: used to extract parameters from indf
  parnames = c("a", "tp", "b", "ga", "p", "la", "q", "er")
  
  #run hitcontinner for every row of indf
  if(mc.cores > 1){
    plan(multiprocess, workers = mc.cores)
    pin = future_lapply(1:nrow(indf), function(i){indf[i, parnames]})
    hitcall = future_mapply(hitcontinner, conc= xs, resp = ys, top = indf$top, cutoff = newcutoff, er = indf$er, ps = pin, 
                     fit_method = indf$fit_method, caikwt = indf$caikwt, mll = indf$mll,
                     future.globals = structure(TRUE, add = c("cnst", "poly1", 
                     "poly2", "pow", "exp2", "exp3", "exp4", "exp5", "hillfn", "gnls") ))
    plan("default")
  } else {
    pin = lapply(1:nrow(indf), function(i){indf[i, parnames]})
    hitcall = mapply(hitcontinner, conc= xs, resp = ys, top = indf$top, cutoff = newcutoff, er = indf$er, ps = pin, 
                     fit_method = indf$fit_method, caikwt = indf$caikwt, mll = indf$mll)
  }

  return(hitcall)

}

#' Continuous Hitcalls Inner
#' 
#' Calculates continuous hitcall using 3 statistical metrics.
#' 
#' This function is called either directly from pathwayConcRespCore_pval or
#' via hitcont using PATHWAY_CR. Details of how to compute function input are in
#' pathwayConcRespCore_pval.
#'
#' @param conc Vector of concentrations.
#' @param resp Vector of responses.
#' @param top Model top.
#' @param cutoff Desired cutoff.
#' @param er Model error parameter.
#' @param ps Vector of used model parameters in order: a, tp, b, ga, p, la, q, er.
#' @param fit_method Name of winning fit method (should never be constant).
#' @param caikwt Aikaike weight of constant model relative to winning model.
#' @param mll Maximum log-likelihood of winning model.
#'
#' @return Continuous hitcall between 0 and 1.
#' @export
#'
#' @examples
#' conc = c(.03,.1,.3,1,3,10,30,100)
#' resp = c(0,.1,0,.2,.6,.9,1.1,1)
#' top = 1.023239
#' er = -3.295307
#' ps = c(1.033239, 2.453014, 1.592714, er = -3.295307) #tp,ga,p,er
#' fit_method = "hill"
#' caikwt = 1.446966e-08
#' mll = 12.71495
#' hitcontinner(conc,resp,top,cutoff = 0.8, er,ps,fit_method, caikwt, mll)
#' hitcontinner(conc,resp,top,cutoff = 1, er,ps,fit_method, caikwt, mll)
#' hitcontinner(conc,resp,top,cutoff = 1.2, er,ps,fit_method, caikwt, mll)
hitcontinner = function(conc, resp, top, cutoff, er, ps, fit_method, caikwt, mll){
  
  #Each P represents the odds of the curve being a hit according to different criteria; multiply all Ps to get hit odds overall
  if(fit_method == "none") return(0)
  if(fit_method == "hill") fname = "hillfn" else fname = fit_method
  
  #caikwt is constant model aikaike weight vs all other models. Represents probability that constant model is correct
  P1 = 1-caikwt
  
  P2 = 1
  for(y in resp){
    #multiply odds of each point falling below cutoff to get odds of all falling below
    P2 = P2*pt((y-sign(top)*cutoff)/exp(er),4, lower.tail = top < 0) #use lower tail for positive top and upper tail for neg top 
  }
  P2 = 1- P2 #odds of at least one point above cutoff
  
  # P3 = pnorm((top-cutoff)/topsd) #odds of top above cutoff
  #assume ps may have nas in them
  ps = ps[!is.na(ps)]
  P3 = toplikelihood(fname, cutoff, conc, resp, ps, top, mll) #odds of top above cutoff
  
  #multiply three probabilities
  return(P1*P2*P3)
  
}
  
#' Top Likelihood
#' 
#' Probability of top being above cutoff.
#' 
#' Should only be called by hitcontinner. Uses profile likelihood, similar
#' to bmdbounds. Here, the y-scale type parameter is substituted in such a
#' way that the top equals the cutoff. Then the log-likelihood is compared to
#' the maximum log-likelihood using chisq function to retrieve probability.
#'
#' @param fname Model function name (equal to model name except hill which
#'   uses "hillfn")
#' @param cutoff Desired cutoff.
#' @param conc Vector of concentrations.
#' @param resp Vector of responses.
#' @param ps Vector of parameters, must be in order: a, tp, b, ga, p, la, q, er
#' @param top Model top.
#' @param mll Winning model maximum log-likelihood.
#'
#' @return Probability of top being above cutoff.
#' @export
#'
#' @examples
#' fname = "hillfn"
#' conc = c(.03,.1,.3,1,3,10,30,100)
#' resp = c(0,.1,0,.2,.6,.9,1.1,1)
#' ps = c(1.033239, 2.453014, 1.592714, er = -3.295307)
#' top = 1.023239
#' mll = 12.71495
#' toplikelihood(fname, cutoff = .8, conc, resp, ps, top, mll)
#' toplikelihood(fname, cutoff = 1, conc, resp, ps, top, mll)
#' toplikelihood(fname, cutoff = 1.2, conc, resp, ps, top, mll)
toplikelihood = function(fname, cutoff, conc, resp, ps, top, mll){
  
  #reparameterize so that top is exactly at cutoff
  if(fname == "exp2"){
    ps[1] = cutoff/( exp(max(conc)/ps[2]) - 1 )
  } else if(fname == "exp3"){
    ps[1] = cutoff/( exp((max(conc)/ps[2])^ps[3]) - 1 )
  } else if(fname == "exp4"){
    ps[1] = cutoff
  } else if(fname == "exp5"){
    ps[1] = cutoff
  } else if(fname == "hillfn"){
    ps[1] = cutoff
  } else if(fname == "gnls"){
    #approximating actual top with theoretical top for convenience.
    ps[1] = cutoff
  } else if(fname == "poly1"){
    ps[1] = cutoff/max(conc)
  } else if(fname == "poly2"){
    ps[1] = cutoff/(max(conc)/ps[2] + (max(conc)/ps[2])^2 ) 
  } else if(fname == "pow"){
    ps[1] = cutoff/(max(conc)^ps[2])
  }
  #get loglikelihood of top exactly at cutoff, use likelihood profile test
  # to calculate probability of being above cutoff
  loglik = tcplObj(p = ps, conc = conc, resp = resp, fname = fname)
  if(abs(top) >= cutoff) out = (1 + pchisq(2*(mll - loglik), 1))/2
  if(abs(top) < cutoff) out = (1 - pchisq(2*(mll - loglik), 1))/2
  
  return(out)
  
}