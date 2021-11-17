#--------------------------------------------------------------------------------------
#
# general_utils.R utilities for managing ToxCast data
#
# August 2017
# Richard Judson
#
# US EPA
# Questions, comments to: judson.richard@epa.gov
#
#
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
printCurrentFunction <- function(comment.string=NA) {
  # prints the name of the current function
  #
  # Args:
  #   comment.string: Adds the value to the string to be printed
  #
  # Returns:
  #   no values returned
  #
  cat("=========================================\n")
	curcall <- sys.call(sys.parent(n=1))[[1]]
	cat(curcall,"\n")
	if(!is.na(comment.string))	cat(comment.string,"\n")
	cat("=========================================\n")
	flush.console()
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
aORb <<- function(V,M) {
  # do an 'or' on a chunk of a matrix
  #
  # Args:
  #   V: first matrix
  #   M: second matrix
  #
  # Returns:
  #   the result of the or operation
  #
  res <- c()
  for(i in 1:dim(M)[1]) res <- c(res,sum(V|M[i,]))
  return(res)
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
colMin <- function(x) {
  # Returns the minimum value for each column of a matrix
  #
  # Args:
  #   x: a numerical matrix or data frame
  #   
  # Returns:
  #   a vector of column minimum values
  #
  ret <- apply(x,FUN=min,MARGIN=2)
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
colMax <- function(x) {
  # Returns the maximum value for each column of a matrix
  #
  # Args:
  #   x: a numerical matrix or data frame
  #   
  # Returns:
  #   a vector of column maximum values
  #
  ret <- apply(x,FUN=max,MARGIN=2)
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
rowMin <- function(x) {
  # Returns the minimum value for each row of a matrix
  #
  # Args:
  #   x: a numerical matrix or data frame
  #   
  # Returns:
  #   a vector of row minimum values
  #
  ret <- apply(x,FUN=min,MARGIN=1)
}
#-----------------------------------------------------------------------------------
#
# minimum by row
#
#-----------------------------------------------------------------------------------
rowMax <- function(x) {
  # Returns the maximum value for each row of a matrix
  #
  # Args:
  #   x: a numerical matrix or data frame
  #   
  # Returns:
  #   a vector of row maximum values
  #
  ret <- apply(x,FUN=max,MARGIN=1)
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#
# median by row
#
#-----------------------------------------------------------------------------------
rowMed <- function(x) {
  # Returns the median value for each row of a matrix
  #
  # Args:
  #   x: a numerical matrix or data frame
  #   
  # Returns:
  #   a vector of row median values
  #
  ret <- apply(x,FUN=median,MARGIN=1)
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
contains <- function(x,query,verbose=F) {
  # detects if a first string contains a second string
  #
  # Args:
  #   x: the first string
  #   query: the second string
  #   verbose: if TRUE, the two strings are printed
  #   
  #
  # Returns:
  #   if x contains qury, retrun TRUE, FALSE otherwise
  #
  if(verbose) {
    print(x)
    print(query)
  }
  if(is.null(x)) return(FALSE)
  if(is.null(query)) return(FALSE)
  if(is.na(x)) return(FALSE)
  if(is.na(query)) return(FALSE)
  x <- stri_trans_tolower(x)
  query <- stri_trans_tolower(query)
  
  val <- sum(grep(query,x,fixed=T))
  if(val>0) return(TRUE)
  else return(FALSE)
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
histLog <- function(x,y,ylim,xlab,ylab,main,cytotox.median, cytotox.min, cytotox.max) {
  # plot a histogram on a log scale Calculate at the hit distribution by chemical
  #
  # Args:
  #   x:
  #   y:
  #   ylim: the limits along the y-axis
  #   xlab: the label for the x-axis
  #   ylab: the label for the y-axis
  #   main: the main label for the plot
  #   cytotx.median: value of the cytotox.median in uM
  #   cytotx.min: value of the cytotox.min in uM
  #   cytotx.max: value of the cytotox.max in uM
  #   
  # Returns:
  #   no values returned
  #
  printCurrentFunction()
  plot(x[1:length(y)],y,type="n",col="gray40",lwd=2.5,log="x",xlab=xlab,ylim=ylim,ylab=ylab,main=main,cex.lab=1.2,cex.axis=1.2)
  if(cytotox.min<100) {
    rect(cytotox.min,ylim[2],max(x),0,col="gray80")
    lines(c(cytotox.median,cytotox.median),ylim,col="red",lwd=3)
  }
  for(i in 1:length(y)) {
    rect(x[i],y[i],x[i+1],0)
  }
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
fixCasrn <- function(casrn,cname="",verbose=F) {
  # Fix a CASRN that has one of several problems
  #
  # Args:
  #   casrn: input CASRN to be fixed
  #   cname: an optional chemical name
  #   verbose: if TRUE, print hte input values
  #
  # Returns:
  #   the fixed CASRN
  #
  if(verbose) cat("input: ",cname,":",casrn,"\n")
  if(contains(casrn,"NOCAS")) return(casrn)
  doit <- T
  while(doit) {
    if(substr(casrn,1,1)=="0") casrn <- substr(casrn,2,nchar(casrn))
    else doit <- F
  }
  
  if(!contains(casrn,"-")) {
    nc <- nchar(casrn)
    ctemp <- casrn
    right <- substr(ctemp,nc,nc)
    mid <- substr(ctemp,nc-2,nc-1)
    left <- substr(ctemp,1,nc-3)
    casrn <- paste(left,"-",mid,"-",right,sep="")
  }
  if(!is.na(cname)) {
    if(cname=="epsilon-Hexachlorocyclohexane (epsilon-HC)") casrn <- "6108-10-7"
    if(cname=="Captafol") casrn <- "2425-06-1"
    if(cname=="Hydrogen sulfide") casrn <- "7783-06-4"
    if(cname=="Picloram") casrn <- "1918-02-1"
    if(cname=="Dodine") casrn <- "2439-10-3"
    if(cname=="Mancozeb") casrn <- "8018-01-7"
  }
  if(verbose) cat("output: ",cname,":",casrn,"\n")
  return(casrn)
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
stringEscape <- function(x) {
  # Escape a string, typically used for chemical names to be put into a database
  #
  # Args:
  #   x: the name to be escaped
  #   
  #
  # Returns:
  #   teh escaped nameno values returned
  #
  x <- str_replace_all(x,"\'","")
  x <- str_replace_all(x,"\"","")
  return(x)
}

