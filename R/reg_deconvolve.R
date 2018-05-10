#' Compute the measurement error version of the Nadaraya-Watson regression
#' estimator
#'
#' Estimates \eqn{m(x) = E[Y | X = x]} from data \eqn{(W, Y)} where 
#' \eqn{W = X + U}. 
#'
#' #' The function \code{reg_deconvolve} chooses from one of two different 
#' methods depending on how the error distribution is defined.
#' 
#' \strong{Error from Replicates:} If both \code{W} and \code{W2} are supplied 
#' then the error is calculated using replicates. This method was prototyped in  
#' Delaigle, Hall, and Meister 2008 and then further refined in Delaigle and  
#' Hall 2016, and Camirand, Carroll, and Delaigle 2018. 
#' 
#' \strong{Homoscedastic Error:} If the errors are defined by either a single 
#' function \code{phiU}, or a single value \code{sd_U} along with its 
#' \code{errortype} then the method used is as described in Fan and Truong 1993.
#' 
#' @param W A vector of the univariate contaminated data W_1, ..., W_n.
#' @param W2 A vector of replicate measurements. If supplied, then the error 
#' will be estimated using replicates.
#' @param Y A vector of the response data Y_1, ..., Y_n.
#' @param xx A vector of x values on which to compute the regression estimator.
#' @param errortype The distribution type of \eqn{U}. Either "laplace" or
#' "normal". If you define the errors this way then you must also provide
#' \code{sd_U} but should not provide \code{phiU}. Argument is case-insensitive
#' and partially matched.
#' @param sd_U The standard deviation of \eqn{U}. This does not need to be
#' provided if you define your error using phiU and provide \code{bw} and
#' \code{rho}.
#' @param phiU A function giving the characteristic function of \eqn{U}. You 
#' should only define the errors this way if you also provide \code{bw} and 
#' \code{rho}. If you define the errors this way then you should not provide 
#' \code{errortype}. 
#' @param bw The bandwidth to use. If you provide this then you should also
#' provide \code{rho}.
#' @param rho The ridge parameter to use. If you provide this then you should
#' also provide \code{bw}.
#' @param n_cores Number of cores to use when calculating the bandwidth. If
#' \code{NULL}, the number of cores to use will be automatically detected.
#' @param kernel_type The deconvolution kernel to use. The default kernel has
#' characteristic function \eqn{(1-t^2)^3}.
#' @param seed Set seed for SIMEX. Otherwise a default seed will be automatically set.
#'
#' @return An object of class deconvolve containing the regression estimator, 
#' as well as the bandwidth and ridge parameter rho. Using SIMEX to choose 
#' smoothing-parameters. See Delaigle and Hall 2008.
#'
#' @section Warnings:
#' \itemize{
#' \item If provided, the bandwidth \code{h} and ridge parameter \code{rho} need
#' to be consistent. You should either provide both or neither.
#' \item The estimator can also be computed using the Fast Fourier Transform,
#' which is faster, but more complex. See Delaigle and Gijbels 2007.
#' }
#'
#' @section References:
#' Fan,  J.,  and Truong,  Y. K. (1993),  Nonparametric Regression With Errors
#' in Variables,  \emph{The Annals of Statistics}.  21,  1900-1925.
#'
#' Delaigle, A. and Gijbels, I. (2007). Frequent problems in calculating
#' integrals and optimizing objective functions: a case study in density
#' deconvolution. \emph{Statistics and Computing}, 17, 349 - 355.
#' 
#' Delaigle, A. and Hall, P. (2008). Using SIMEX for smoothing-parameter choice
#' in errors-in-variables problems. \emph{Journal of the American Statistical
#' Association}, 103, 481, 280-287
#' 
#' Delaigle, A., Hall, P., and Meister, A. (2008). On Deconvolution with  
#' repeated measurements. \emph{Annals of Statistics}, 36, 665-685 
#' 
#' Delaigle, A. and Hall, P. (2016). Methodology for non-parametric 
#' deconvolution when the error distribution is unknown. \emph{Journal of the 
#' Royal Statistical Society: Series B (Statistical Methodology)}, 78, 1, 
#' 231-252.
#' 
#' Camirand, F., Carroll, R.J., and Delaigle, A. (2018). Estimating the  
#' distribution of episodically consumed food measured with errors.  
#' \emph{Manuscript.} 
#' 
#' @author Aurore Delaigle, Timothy Hyndman, Tianying Wang
#' 
#' @example man/examples/reg_deconvolve_eg.R
#'
#' @export

reg_deconvolve <- function(W, Y, W2 = NULL, xx = seq(min(W), max(W), length.out = 100), 
                           errortype = NULL, sd_U = NULL, phiU = NULL, 
                           bw = NULL, rho = NULL, n_cores = NULL, 
                           kernel_type = c("default", "normal", "sinc"), 
                           seed = NULL) {

    # Partial matching ---------------------------------------------------------
    dist_types <- c("normal", "laplace")
    if (!is.null(errortype)) {
        errortype <- dist_types[pmatch(tolower(errortype), dist_types)]
        if (is.na(errortype)) {
            stop("Please provide a valid errortype.")
        }
    }

    kernel_type <- match.arg(kernel_type)

    if (!is.null(W2)) {
        errors <- "rep"
    } else {
        errors <- "hom"
    }

    # Check inputs -------------------------------------------------------------
    if (is.null(errortype) & is.null(phiU) & is.null(W2)) {
        stop("You must provide either errortype, phiU, or W2.")
    }
    
    if (errors == "hom") {
        if (!is.null(errortype) & is.null(sd_U)) {
            stop("You must provide sd_U along with errortype.")
        }

        if ((is.null(bw) | is.null(rho)) & (is.null(sd_U) | is.null(errortype))) {
            stop("If the bandwidth is not provided then you must provide errortype and sd_U.")
        }

        if ((is.null(bw) | is.null(rho)) & is.null(sd_U)){
            stop("You must provide sd_U if you do not provide bw and rho.")
        }
    }

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

    # Convert errortype to phiU ------------------------------------------------
    if (is.null(phiU)){
        phiU <- create_phiU(errors = "hom", errortype, sd_U)
    }

    # Calculate bandwidth ------------------------------------------------------
    if (is.null(bw) | is.null(rho)) {
        outcome_tmp <- bandwidth(W = W, 
                                 W2 = W2,
                                 errortype = errortype, 
                                 sd_U = sd_U,
                                 phiU = phiU, 
                                 Y = Y, 
                                 algorithm = "SIMEX", 
                                 n_cores = n_cores, 
                                 kernel_type = kernel_type, 
                                 seed = seed)
        bw <- outcome_tmp$h
        rho <- outcome_tmp$rho
    }

    # Compute estimate for m(X) ------------------------------------------------
    y <- NWDecUknown(xx, W, Y, phiU, bw, rho, phiK, tt, deltat)

    structure(list(pdf = y, bw = bw, rho = rho, x = xx), 
                   class = c("reg_deconvolve", "list"))


}
