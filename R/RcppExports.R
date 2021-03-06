# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title The probability density of standard Laplace distribution using Rcpp
#' @description The probability density of standard Laplace distribution using Rcpp
#' @param x the value of the variable
#' @return the probability density of the value of x
#' @examples
#' \dontrun{
#' f(4)
#' }
#' @export
f <- function(x) {
    .Call('_StatComp20001_f', PACKAGE = 'StatComp20001', x)
}

#' @title A generation of standard Laplace distribution using Rcpp
#' @description A generation of standard Laplace distribution using Rcpp
#' @param sigma the parameter of proposal distribution
#' @param x0 the initial values
#' @param N the number of samples
#' @return a random sample of size \code{N}
#' @examples
#' \dontrun{
#' LaplaceC <- rwMetropolis(0.2,5,3000)
#' }
#' @export
rwMetropolis <- function(sigma, x0, N) {
    .Call('_StatComp20001_rwMetropolis', PACKAGE = 'StatComp20001', sigma, x0, N)
}

