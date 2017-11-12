#' Compute the measurement error version of the Nadaraya-Watson regression
#' estimator
#'
#' Goal: estimate m where Y=m(X)+epsilon,  and we observe data on (W, Y),  where
#' W=X+U.
#' See Fan,  J.,  and Truong,  Y. K. (1993),  Nonparametric Regression With
#' Errors in Variables,  The Annals of Statistics,  21,  1900-1925

#' @param xx vector of x-values where to compute the regression estimator
#' @param W vector of contaminated data W_1, ..., W_n
#' @param Y vector of data Y_1, ..., Y_n
#' in errors-in-variables problems.  JASA,  103,  280-287.
#' @param errortype 'Lap' for Laplace errors and 'norm' for normal errors.
#' @param sigU parameter of Laplace or normal errors used only to define characteristic function of the error.
#' @param n_cores cores used to do parallel computation to calculate bandwidth.
#'
#' @return Regression estimator, bandwidth and ridge parameter rho. See Delaigle,  A. and Hall,  P. (2008). Using SIMEX for smoothing-parameter choice
#'
#' @section Warnings:
#' \itemize{
#' \item h and rho need to be consistent. They need to be both specified or unspecified.
#' \item The range of t-values -1 and 1 correspond to the support of phiK.
#' \item If the grid of t-values is fine enough,  the estimator can simply be
#' computed like here without having problems with oscillations.
#' \item The phiK here is also used to compute the bandwidth.
#' \item The estimator can also be computed using the Fast Fourier Transform,  which
#' is faster,  but more complex.
#' \item See Delaigle,  A. and Gijbels,  I. (2007). Frequent problems in calculating
#' integrals and optimizing objective functions: a case study in density
#' deconvolution.   Statistics and Computing,   17,   349 - 355.}
#'
#' @section References:
#' Fan,  J.,  and Truong,  Y. K. (1993),  Nonparametric Regression With Errors in Variables,  \emph{The
#' Annals of Statistics}.  21,  1900-1925.
#'
#' Delaigle,  A. and Hall,  P. (2008). Using SIMEX for smoothing-parameter choice in errors-in-variables
#' problems.  \emph{JASA},  103,  280-287.
#'
#' @example man/examples/reg_deconvolve_eg.R
#'
#' @export

reg_deconvolve <- function(xx, W, Y, errortype, sigU, n_cores = NULL) {

    # --------------------------------------------------------
    # Preliminary calculations and initialisation of functions
    # --------------------------------------------------------
    W <- as.vector(W)
    n <- length(W)
    # error can only be normal or laplace
    if (errortype == "Lap") {
        phiU <- function(t) {
            1 / (1 + sigU^2 * t^2)
        }
    }
    if (errortype == "norm") {
        phiU <- function(t) {
            exp(-sigU^2 * t^2 / 2)
        }
    }

    # phiK: Fourier transform of the kernel K. You can change this if you wish
    # to use another kernel but make sure
    # you change the range of t-values,  which should correspond to the support
    # of phiK

    phiK <- function(t) {
        (1 - t^2)^3
    }


    # Range of t-values (must correspond to the domain of phiK)
    dt <- .0002
    tt <- seq(-1, 1, dt)
    longt <- length(tt)
    dim(tt) <- c(length(tt), 1)

    outcome_tmp = deconvolve::bandwidth(W = W, Y = Y, errortype = errortype, sigU = sigU, algorithm = "SIMEX", n_cores = n_cores)
        h = outcome_tmp$h
        rho = outcome_tmp$rho
    # Compute the empirical characteristic function of W (times n) at t/h:
    # \hat\phi_W(t/h)
    OO <- t(outerop(tt/h, t(W),"*"))
    csO <- cos(OO)
    snO <- sin(OO)
    rm(OO)

    rehatphiW <- apply(csO, 2, sum)
    imhatphiW <- apply(snO, 2, sum)


    # Compute \sum_j Y_j e^{itW_j/h}
    dim(Y) <- c(1, n)
    renum <- Y %*% csO
    imnum <- Y %*% snO


    # Compute numerator and denominator of the estimator separately
    # Numerator: real part of
    # (2*pi*h)^(-1) \int e^{-itx/h} n^(-1)\sum_j Y_j e^{itW_j/h} \phi_K(t)/\phi_U(t/h) dt
    #
    # Denominator: real part of
    # (2*pi*h)^(-1) \int e^{-itx/h} \hat\phi_W(t/h) \phi_K(t)/\phi_U(t/h) dt


    xt <- outerop(tt / h, t(xx),"*")
    cxt <- cos(xt)
    sxt <- sin(xt)
    rm(xt)

    phiUth <- phiU(tt / h)
    matphiKU <- phiK(tt) / phiUth
    dim(matphiKU) <- c(1, longt)

    Den <- (rehatphiW * matphiKU) %*% cxt + (imhatphiW * matphiKU) %*% sxt
    Num <- (renum * matphiKU) %*% cxt + (imnum * matphiKU) %*% sxt
    Num <- Num * (dt / (n * h * 2 * pi))
    Den <- Den * (dt / (n * h * 2 * pi))


    # If denomintor is too small,  replace it by the ridge rho
    dd <- Den
    dd[which(dd < rho)] <- rho

    # Finally obtain the regression estimator
    y <- Num / dd

    return(list(y=as.vector(y), bw=h, rho = rho))
}



#OUTEROP calculate an outer operation on two vectors
#   function y=outerop(A,B,OPERATOR)
#
#   Calculates resultant matrix when the OPERATOR is applied
#   to all combinations of the elements of vector A and the
#   elements of vector B e.g. the outer product of A and B
#   is outerop(A,B,'*'), the outer sum of A and B
#   is outerop(A,B,'+')
#
#   If OPERATOR is omitted '+' is assumed
#
#   This function is equivalent to the
#   APL language's circle.dot operation. e.g.
#   in APL Ao.*B is the outer product
#   of A and B % and Ao.+B is the outer sum.
#   Ideally it would work on matrices but in practice
#   I've only ever needed to use it with vectors.
#
#   Copyright Murphy O'Brien 2005
#   all rights unreserved
#
outerop<-function(a,b,operator)
{
    if (missing(operator))
    {operator="+"}                      # for only two arguments assume outerproduct

    if (operator=="*")                      # common operator
    {y=a%*%b} else
    {outera=as.matrix(a)%*%rep(1,length(b))          # these are the matrices that
    outerb=as.matrix(rep(1,length(a)))%*%b # meshgrid(A,B) would generate
    functionHandle=match.fun(operator)    # new R14 functionality
    y=functionHandle(outera,outerb)  }    # allows faster/neater method
    return (y)
}
