#' #' options set upon package loading
#' .onLoad <- function(lib, pkg) {
#'
#'   utils::globalVariables(c("FCMAT1","FCMAT2_lowconc","GENELISTS","CATALOG","COVMAT","SIGSCOREMAT","DSSTOX"))
#'   utils::globalVariables(c("MAT", "FCMAT1.1", "FCMAT1.2", "CHEM_DICT", "FCMAT2", "GMAT", "SIG_CR", "TOXCAST", "l2fc"))
#'   utils::globalVariables(c("gene","conc","sample_id","dtxsid","casrn","name","time","bmed","cutoff","onesd","proper_name"))
#'   utils::globalVariables(c("bmr","er","fit_method","ac50","bmd","top","top_over_cutoff","hitcall","bmdl","bmdu","conc"))
#'   utils::globalVariables(c("signature_score","size","nsig","super_target","signature","signature_size"))
#'
#'
#'   FCMAT1 <<- NULL
#'   FCMAT2_lowconc <<- NULL
#'   GENELISTS <<- NULL
#'   CATALOG <<- NULL
#'   COVMAT <<- NULL
#'   SIGSCOREMAT <<- NULL
#'   DSSTOX <<- NULL
#'   MAT <<- NULL
#'   FCMAT1.1 <<- NULL
#'   FCMAT1.2 <<- NULL
#'   CHEM_DICT <<- NULL
#'   FCMAT2 <<- NULL
#'   GMAT <<- NULL
#'   SIG_CR <<- NULL
#'   TOXCAST <<- NULL
#'   MAT2 <<- NULL
#'   RES <<- NULL
#' }

