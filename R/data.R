#' Chemical Dictionary Example
#'
#' Chemical dictionary for 10 randomly chosen chemicals from the phase 1 screen.
#'
#' @format A data frame with 80 rows and 7 variables: sample_key (sample id + conc),
#'   sample_id (unique identifier for each sample), conc (concentration), time
#'   (length of the experiment in hours), casrn, name (chemical name),
#'   dtxsid.
#"CHEM_DICT"

#' Endocrine Pathways Example
#'
#' Pathway dictionary for 44 ER/AR related pathways
#'
#' @format A data frame with 44 rows and 8 variables: pathset, pathway, ngene,
#' ngene_in_httr, pathway_class, pathway_super_class, url, gene_list.
#"ENDOCRINE_PATHWAYS"

#' Fold Change Matrix Example
#'
#' Fold change matrix for 10 randomly selected chemicals from the phase 1 screen.
#'
#' @format A sample by gene matrix with 80 rows and 10,341 columns.
#"FCMAT2"

#' Endocrine Pathway Data
#'
#' Pathway data example for 44 ER/AR related pathways
#'
#' @format A named list containing 44 elements. Each elements corresponds to a
#'   pathway and consists of a vector of gene names.
#"pathway_data"
