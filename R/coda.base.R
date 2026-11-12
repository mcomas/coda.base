#' coda.base
#'
#' A minimum set of functions to perform compositional data analysis
#' using the log-ratio approach introduced by John Aitchison (1982)
#' \doi{10.1111/j.2517-6161.1982.tb01195.x}. Main functions
#' have been implemented in C++ for better performance.
#'
#' @docType package
#' @author Marc Comas-Cufí
#' @import Rcpp stats
#' @useDynLib coda.base, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom grDevices dev.flush dev.hold
#' @importFrom graphics par points polygon rect segments strheight strwidth text
#'
#' @name coda.base
"_PACKAGE"
