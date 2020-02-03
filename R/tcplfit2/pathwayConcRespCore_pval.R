#' Pathway Concentration Response Core
#'
#' Core of concentration response curve fitting for pvalue based cutoff. This
#' function calls httrFit to get curve fits, chooses the winning model, extracts
#' the top and ac50, computes the hitcall, and calculates bmd/bmdl/bmdu among other
#' statistics. Nested model selection is used to choose between poly1/poly2, then
#' the model with the lowest AIC (or AICc) is declared the winner. Continuous
#' hitcalls requires httrFit to be run with force.fit = T and "cnst" never to
#' be chosen as the winner.
#'
#' @param row A named list that must include:
#'   \itemize{
#'     \item conc - list of concentrations (not in log units)
#'     \item resp - list of corresponding responses
#'     \item bmed - median of noise estimate.
#'     \item cutoff - noise cutoff
#'     \item onesd - 1 standard deviation of the noise (for bmd calculation)
#'   }
#'   Other elements (usually identifiers, like casrn) of row will be attached to
#'   the final output.
#' @param fitmodels Vector of model names to use.
#' @param conthits conthits = T uses continuous hitcalls, otherwise they're
#'   discrete.
#' @param aicc aicc = T uses corrected AIC to choose winning method; otherwise
#'   regular AIC.
#'
#' @return One row dataframe containing all CR output and statistics and any
#'   identifiers from row.
#' @export
#'
#' @examples
#' conc = list(.03,.1,.3,1,3,10,30,100)
#' resp = list(0,.2,.1,.4,.7,.9,.6, 1.2)
#' row = list(conc = conc, resp = resp, bmed = 0, cutoff = 1, onesd = .5)
#' pathwayConcRespCore_pval(row, conthits = TRUE)
#' pathwayConcRespCore_pval(row, aicc = TRUE)
pathwayConcRespCore_pval <- function(row,
                                     fitmodels = c("cnst", "hill", "gnls", "poly1", "poly2", "pow", "exp2", "exp3",
                                                   "exp4", "exp5"),
                                     conthits = T, aicc = F) {

  #row needs to include cutoff and bmed
  #unpack row into the local environment, for ease: sample_id, dtxsid, casrn, name, time, pathway, size, con, resp
  list2env(row,envir=environment())
  resp = unlist(resp)
  conc = unlist(conc)

  #prepare input
  resp = resp - bmed
  conc = conc[!is.na(resp)]
  resp = resp[!is.na(resp)]

  #run the fits
  if(conthits) fitmodels = unique(c("cnst", fitmodels)) #cnst models must be present for conthits but not chosen
  params <- httrFit(conc, resp, cutoff, force.fit = conthits, bidirectional = T, fitmodels = fitmodels)

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
  conc <- paste(conc,collapse="|")
  resp <- paste(resp,collapse="|")

  #PATHWAY_CR contains the specified columns and any identifying, unused columns
  #that were in pathscoremat/
  identifiers = row[!names(row) %in% c("conc", "resp", "bmed", "onesd", "cutoff")]
  name.list <- c("n_gt_cutoff","cutoff", "fit_method",
                 "top_over_cutoff", "rmse", "a", "b", "tp", "p", "q", "ga", "la", "er", "bmr", "bmdl", "bmdu", "caikwt",
                 "mll","hitcall", "ac50","ac50_loss","top", "ac5","ac10","ac20", "acc", "ac1sd", "bmd", "conc", "resp")
   row = as.data.frame(c(identifiers, mget(name.list)), stringsAsFactors = F)
  return(row)
}


#' Nest Select
#'
#' Chooses between nested models.
#'
#' @param aics Named vector of model aics (can include extra models).
#' @param mod1 Name of model 1, the model with fewer degrees of freedom.
#' @param mod2 Name of model 2, the model with more degrees of freedom.
#' @param dfdiff Absolute difference in number of degrees of freedom
#'   (i.e. the difference in parameters).
#' @param pval P-value for nested model test.
#'
#' @return Named aic vector with losing model removed.
#' @export
#'
#' @examples
#' aics = c(-5,-6,-3)
#' names(aics) = c("poly1", "poly2", "hill")
#' nestselect(aics, "poly1", "poly2", 1)
#'
#' aics = c(-5,-7,-3)
#' names(aics) = c("poly1", "poly2", "hill")
#' nestselect(aics, "poly1", "poly2", 1)
nestselect = function(aics, mod1, mod2, dfdiff, pval = .05){


  if(is.na(aics[mod1])){
    loser = mod1 #if model 1 AIC is NA it is the loser
  } else if(isTRUE(aics[mod2] <= aics[mod1])){
    #if both aics exist and model 2 aic is lower and model 2 passes the ratio test,
    #model 2 wins and model 1 loses; otherwise, model 2 loses.
    chisq = aics[mod1] - aics[mod2] + 2*dfdiff #2*loglikelihood(poly2) - 2*loglikelihood(poly1)
    ptest =  1-pchisq(chisq, dfdiff)
    if(ptest < pval) loser = mod1 else loser = mod2
  } else loser = mod2

  return(aics[names(aics) != loser])
}
