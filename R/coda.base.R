#' coda.base
#'
#' A minimum set of functions to perform compositional data analysis
#' using the log-ratio approach introduced by John Aitchison (1982)
#' <http://www.jstor.org/stable/2345821>. Main functions
#' have been implemented in c++ for better performance.
#'
#' @docType package
#' @author Marc Comas-Cufí
#' @import Rcpp stats
#' @useDynLib coda.base, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @name coda.base
NULL
