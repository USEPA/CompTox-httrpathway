#' BMD Bounds
#' 
#' Computes BMDU or BMDL.
#' 
#' Takes in concentration response fit details and outputs a bmdu or bmdl, as
#' desired. If bmd is not finite, returns NA. If the objective function doesn't
#' change sign or the root finding otherwise fails, it returns NA. These 
#' failures are not uncommon since some curves just don't reach the desired
#' confidence level.
#'
#' @param fit_method Fit method: "exp2", "exp3", "exp4", "exp5", "hill", "gnls",
#'   "poly1", "poly2", or "pow".
#' @param bmr Benchmark response.
#' @param pars Named vector of model parameters: a,b,tp,ga,p,la,q,er output by
#'   httrfit, and in that order.
#' @param conc Vector of concentrations (NOT in log units).
#' @param resp Vector of responses corresponding to given concentrations.
#' @param onesidedp The one-sided p-value. Default of .05 corresponds to 5 
#'   percentile BMDL, 95 percentile BMDU, and 90 percent CI.
#' @param bmd Can optionally input the bmd when already known to avoid 
#'   unnecessary calculation.
#' @param which.bound Returns BMDU if which.bound = "upper"; returns BMDL if 
#'   which.bound = "lower".
#'
#' @return Returns either the BMDU or BMDL.
#' @export
#'
#' @examples
#' conc = c(.03, .1, .3, 1, 3, 10, 30, 100)
#' resp = c(.1,-.1,0,1.1,1.9,2,2.1,1.9)
#' pars = c(tp = 1.973356, ga = 0.9401224, p = 3.589397, er =  -2.698579)
#' bmdbounds(fit_method = "hill", bmr = .5, pars, conc, resp)
#' bmdbounds(fit_method = "hill", bmr = .5, pars, conc, resp, which.bound = "upper")

bmdbounds = function(fit_method, bmr, pars, conc, resp, onesidedp = .05, bmd = NULL, which.bound = "lower"){
  
  #calculate bmd, if necessary
  if(is.null(bmd)) bmd = acy(bmr, as.list(pars), type = fit_method)
  if(!is.finite(bmd)) return(NA_real_)
  
  # hill model's function name is hillfn, other models are not changed.
  if(fit_method == "hill") fname = paste0(fit_method, "fn") else fname = fit_method
  maxloglik = tcplObj(p = pars, conc = conc, resp = resp, fname = fname)
  
  #search for bounds to ensure sign change
  if(which.bound == "lower") {
    xs = 10^seq(-5,log10(bmd), length.out = 100)
    ys = sapply(X = xs, FUN = bmdobj, fname = fname, bmr = bmr, conc = conc, resp = resp, ps = pars, mll = maxloglik,
                onesp = onesidedp, partype = 2)
    if(!any(ys >= 0, na.rm = T) | !any(ys < 0, na.rm = T)) return(NA_real_)
    bmdrange = c(max(xs[ys >= 0]), bmd)
  }
  if(which.bound == "upper") {
    if(fit_method == "gnls"){
      toploc = acy(bmr, as.list(pars), type = "gnls", returntoploc = T)
      xs = 10^seq(log10(bmd), log10(toploc), length.out = 100)
    } else xs = 10^seq(log10(bmd), 5, length.out = 100)
    ys = sapply(X = xs, FUN = bmdobj, fname = fname, bmr = bmr, conc = conc, resp = resp, ps = pars, mll = maxloglik,
                onesp = onesidedp, partype = 2)
    if(!any(ys >= 0, na.rm = T) | !any(ys < 0, na.rm = T)) return(NA_real_)

    bmdrange = c(bmd, min(xs[ys >= 0]))
  }
  
  #use type 2 param. only
  out = try(uniroot(bmdobj, bmdrange, fname = fname, bmr = bmr, conc = conc, resp = resp, ps = pars, mll = maxloglik,
                onesp = onesidedp, partype = 2)$root)
  
  #sometimes there's no lower/upper bound because the model is such a poor fit, in this case, return NA
  if(class(out) == "try-error") return(NA_real_) else return(out)

}

#' BMD Objective Function
#' 
#' Utility function for bmdbounds
#'
#' @param bmd Benchmark dose.
#' @param fname Function name: "exp2", "exp3", "exp4", "exp5", "hillfn", "gnls",
#'   "poly1", "poly2", or "pow".
#' @param bmr Benchmark response.
#' @param conc Vector of concentrations NOT in log units.
#' @param resp Vector of corresponding responses.
#' @param ps Named list of paramters.
#' @param mll Maximum log-likelihood of winning model.
#' @param onesp One-sided p-value.
#' @param partype Number for parameter type. Type 1 is y-scaling: a or tp.
#'   Type 2 is x-scaling: b or ga, when available, a otherwise. Type 3 is
#'   power scaling: p when available, then b or ga, then a if no others.
#'   Since bmd is linked to the x-scale, type 2 should always be used. Other
#'   types can also be vulnerable to underflow/overflow.
#'
#' @return Objective function value to find the zero of.
#' @export
bmdobj= function(bmd, fname, bmr, conc, resp, ps, mll, onesp, partype = 2){
  
  #implements the BMD substitutions in Appendix A of the Technical Report.
  #Changes one of the existing parameters to an explicit bmd parameter through
  # the magic of algebra.
  if(fname == "exp2"){
    if(partype == 1) ps["a"] = bmr/( exp(bmd/ps["b"]) - 1 )
    if(partype == 2) ps["b"] = bmd/( log(bmr/ps["a"] + 1) )
    if(partype == 3) ps["b"] = bmd/( log(bmr/ps["a"] + 1) )
  } else if(fname == "exp3"){
    if(partype == 1) ps["a"] = bmr/( exp((bmd/ps["b"])^ps["p"]) - 1 )
    if(partype == 2) ps["b"] = bmd/( log(bmr/ps["a"] + 1) )^(1/ps["p"])
    if(partype == 3) ps["p"] = log(log(bmr/ps["a"] + 1))/log(bmd/ps["b"])
  } else if(fname == "exp4"){
    if(partype == 1) ps["tp"] = bmr/( 1-2^(-bmd/ps["ga"]))
    if(partype == 2) ps["ga"] = bmd/( -log2(1-bmr/ps["tp"]) )
    if(partype == 3) ps["ga"] = bmd/( -log2(1-bmr/ps["tp"]) )
  } else if(fname == "exp5"){
    if(partype == 1) ps["tp"] = bmr/( 1-2^(-(bmd/ps["ga"])^ps["p"] ))
    if(partype == 2) ps["ga"] = bmd/(( -log2(1-bmr/ps["tp"]) )^(1/ps["p"]))
    if(partype == 3) ps["p"] = log( -log2( 1 - bmr/ps["tp"]) )/log(bmd/ps["ga"])
  } else if(fname == "hillfn"){
    if(partype == 1) ps["tp"] = bmr*( 1 + (ps["ga"]/bmd)^ps["p"])
    if(partype == 2) ps["ga"] = bmd* (ps["tp"]/bmr - 1)^(1/ps["p"])
    if(partype == 3) ps["p"] = log(ps["tp"]/bmr - 1)/log(ps["ga"]/bmd)
  } else if(fname == "gnls"){
    #gnls bounds don't include loss part parameterization
    if(partype == 1) ps["tp"] = bmr * ( (1 + (ps["ga"]/bmd)^ps["p"]) * ( 1 + (bmd/ps["la"])^ps["q"]) )
    if(partype == 2) ps["ga"] = bmd * ( ps["tp"]/(bmr*(1 + (bmd/ps["la"])^ps["q"])) - 1)^(1/ps["p"])
    if(partype == 3) ps["p"] = log(ps["tp"]/(bmr*(1 + (bmd/ps["la"])^ps["q"])) - 1) / log(ps["ga"]/bmd)
  } else if(fname == "poly1"){
    if(partype == 1) ps["a"] = bmr/bmd
    if(partype == 2) ps["a"] = bmr/bmd
    if(partype == 3) ps["a"] = bmr/bmd
  } else if(fname == "poly2"){
    if(partype == 1) ps["a"] = bmr/(bmd/ps["b"] + (bmd/ps["b"])^2 ) 
    if(partype == 2) ps["b"] = 2*bmd/(sqrt(1 + 4*bmr/ps["a"]) - 1)
    if(partype == 3) ps["b"] = 2*bmd/(sqrt(1 + 4*bmr/ps["a"]) - 1)
  } else if(fname == "pow"){
    if(partype == 1) ps["a"] = bmr/(bmd^ps["p"])
    if(partype == 2) ps["a"] = bmr/(bmd^ps["p"])
    if(partype == 3) ps["p"] = log(bmr/ps["a"])/log(bmd)
  }
  
  loglik = tcplObj(p = ps, conc = conc, resp = resp, fname = fname)
  #for bmd bounds, we want the difference between the max log-likelihood and the
  #bounds log-likelihood to be equal to chi-squared at 1-2*onesp (typically .9)
  #with one degree of freedom divided by two. 
  return(mll - loglik - qchisq(1-2*onesp,1)/2)
  
}

