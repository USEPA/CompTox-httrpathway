library(RMySQL)
library(DBI)
#--------------------------------------------------------------------------------------
#' set SQL connection to the database
#' @param server SQL server on which relevant database lives
#' @param user SQL username to access database
#' @param password SQL password corresponding to username
#--------------------------------------------------------------------------------------
setDBConn <- function(server="mysql-res1.epa.gov",user="rjudson",password=NA) {
  printCurrentFunction()
  DB.SERVER <<- server
  DB.USER <<- user
  DB.PASSWORD <<- password
}  











