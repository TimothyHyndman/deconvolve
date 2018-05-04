#' Compute the measurement error version of the Nadaraya-Watson regression
#' estimator when the errors are known.
#'
#' Estimates \eqn{m(x) = E[Y | X = x]} from data \eqn{(W, Y)} where 
#' \eqn{W = X + U} and the \eqn{U} are homoscedastic and known. See Fan and 
#' Truong 1993.
#'
#' @param W A vector of the univariate contaminated data W_1, ..., W_n.
#' @param Y A vector of the response data Y_1, ..., Y_n.
#' @param xx A vector of x values on which to compute the regression estimator.
#' @param phiU A function giving the characteristic function of \eqn{U}.
#' @param h The bandwidth to use.
#' @param rho The ridge parameter to use.
#' @param kernel_type The deconvolution kernel to use. The default kernel has
#' characteristic function \eqn{(1-t^2)^3}.
#'
#' @return An object of class reg_deconvolve containing the regression estimator, 
#' as well as the bandwidth and ridge parameter rho. Using SIMEX to choose 
#' smoothing-parameters. See Delaigle and Hall 2008.
#'
#' @section Warnings:
#' \itemize{
#' \item The bandwidth \code{h} and ridge parameter \code{rho} need to be 
#' consistent. You can calculate these using the function \code{bandwidth}.
#' \item The estimator can also be computed using the Fast Fourier Transform,
#' which is faster, but more complex. See Delaigle and Gijbels 2007.
#' }
#'
#' @section References:
#' Fan,  J.,  and Truong,  Y. K. (1993),  Nonparametric Regression With Errors
#' in Variables,  \emph{The Annals of Statistics}.  21,  1900-1925.
#'
#' Delaigle, A. and Hall, P. (2008). Using SIMEX for smoothing-parameter choice
#' in errors-in-variables problems. \emph{Journal of the American Statistical
#' Association}, 103, 481, 280-287
#'
#' Delaigle, A. and Gijbels, I. (2007). Frequent problems in calculating
#' integrals and optimizing objective functions: a case study in density
#' deconvolution. \emph{Statistics and Computing}, 17, 349 - 355.
#' 
#' @author Aurore Delaigle, Timothy Hyndman, Tianying Wang
#' 
#' @example man/examples/hom_regression_U_known_eg.R
#'
#' @export

hom_regression_U_known <- function(W, 
                                   Y, 
                                   phiU, 
                                   h, 
                                   rho,
                                   xx = seq(min(W), max(W), length.out = 100),  
                                   kernel_type = c("default", "normal", "sinc")) {

    kernel_type <- match.arg(kernel_type)

    # Check inputs -------------------------------------------------------------
    if (kernel_type == "normal") {
        warning("You should only use the 'normal' kernel when the errors are 
            Laplace or convolutions of Laplace.")
    }

    if (kernel_type == "sinc") {
        warning("You should ensure that you are not using a plug-in bandwidth 
            method for the bandwidth.")
    }

    # --------------------------------------------------------------------------
    kernel_list <- kernel(kernel_type)
    phiK <- kernel_list$phik
    tt <- kernel_list$tt
    deltat <- tt[2] - tt[1]

    # Compute estimate for m(X) ------------------------------------------------
    y <- NWDecUknown(xx, W, Y, phiU, h, rho, phiK, tt, deltat)

    structure(list(pdf = y, h = h, rho = rho, x = xx), 
                   class = c("reg_deconvolve", "list"))

}
