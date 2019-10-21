#' R Squared
#' 
#' Calculate coefficient of determination.
#' 
#' Note that order matters: R2(x,y) != R2(y,x) in general.
#'
#' @param y Vector of actual values. 
#' @param pred Vector of corresponding predicted values.
#'
#' @return Coefficient of determination.
#' @export
#'
#' @examples
#' R2(c(1:10), c(1:10*.8))
#' R2(c(1:10*.8), c(1:10))
R2 = function(y,pred){
  num = sum((y - pred)^2, na.rm = T) 
  denom = sum((y-mean(y, na.rm = T))^2,na.rm = T)
  return(1 - num/denom)
}

#' Weighted Root-mean-square-error
#' 
#' Computes root-mean-square error with weighted average.
#' 
#' x,y,w should all be the same length. Order of x and y won't change output.
#' 
#'
#' @param x First vector of numbers.
#' @param y Second vector of numbers.
#' @param w Vector of weights.
#'
#' @return Weighted RMSE.
#' @export
#'
#' @examples
#' WRMSE(1:3, c(1,3,5), 1:3)
WRMSE = function(x,y,w){
  return(sqrt(sum(w*(x-y)^2,na.rm = T)/sum(w)))
  
}

#' Root-mean-square-error
#' 
#' Computes root-mean-square-error between two vectors.
#'
#' @param x First vector.
#' @param y Second vector.
#'
#' @return RMSE
#' @export
#'
#' @examples
#' RMSE(1:3, c(1,3,5))
RMSE = function(x,y){
  return(sqrt(mean((x-y)^2,na.rm = T)))
  
}

#' Area Under the Curve
#' 
#' Compute AUC for an ROC curve.
#' 
#' Uses trapezoid rule numerical integration to approximate AUC. Will be more
#' accurate with more fine-grained inputs.
#'
#' @param tpr Vector of true positive rates.
#' @param fpr Vector of false positive rates.
#'
#' @return AUC
#' @export
#'
#' @examples
#' auc(c(0,.5,1), c(0,.5,1))
#' auc(c(0,1,1), c(0,.5,1))
auc = function(tpr, fpr){
  #trapezoid rule integration
  auc = 0
  for(i in 1:(length(tpr)-1)){
    auc = auc + (tpr[i] + tpr[i+1])*(fpr[i+1]-fpr[i])/2
  }
  return(abs(auc))
}