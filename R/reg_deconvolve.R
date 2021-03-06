#' Compute the measurement error version of the Nadaraya-Watson regression
#' estimator
#'
#' Estimates \eqn{m(x) = E[Y | X = x]} from data \eqn{(W, Y)} where 
#' \eqn{W = X + U}. 
#'
#' #' The function \code{reg_deconvolve} chooses from one of two different 
#' methods depending on how the error distribution is defined.
#' 
#' \strong{Error from Replicates:} If both \code{W1} and \code{W2} are supplied 
#' then the error is calculated using replicates. This method was prototyped in  
#' Delaigle, Hall, and Meister 2008 and then further refined in Delaigle and  
#' Hall 2016, and Camirand, Carroll, and Delaigle 2018. 
#' 
#' \strong{Homoscedastic Error:} If the errors are defined by either a single 
#' function \code{phiU}, or a single value \code{sd_U} along with its 
#' \code{errortype} then the method used is as described in Fan and Truong 1993.
#' 
#' The order in which we choose the methods is as follows:
#' \enumerate{
#'  \item If provided, use \code{phiU} to define the errors, otherwise
#'  \item If provided use \code{errortype} and \code{sd_u} to define the errors, otherwise
#'  \item If provided, use the vector of replicates \code{W2} to estimate the error distribution.
#' }
#' 
#' Note that in both 1 and 2, if a vector of replicates \code{W2} is provided we
#' augment the data in \code{W1} with that in \code{W2}.
#' 
#' @param Y A vector of the response data Y_1, ..., Y_n.
#' @param xx A vector of x values on which to compute the regression estimator.
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
#' @param kernel_type The deconvolution kernel to use. The default kernel has characteristic function 
#' \eqn{(1-t^2)^3} for \eqn{t \in [-1,1]}. The normal kernel is the standard normal density.
#' The sinc kernel has characteristic function equal to 1 for \eqn{t \in [-1,1]}
#' @param seed Set seed for SIMEX. Allows for reproducible results using SIMEX.
#' Otherwise a default seed will be automatically set.
#' @param use_alt_SIMEX_rep_opt Only used with SIMEX using replicates. If 
#' \code{TRUE}, performs SIMEX on \eqn{W = (W1 + W2)/2} and samples \eqn{U*} 
#' from (W1 - W2). The default performs SIMEX on \eqn{W = (W1, W2)} and 
#' and samples \eqn{U*} from \eqn{(W1 - W2)/\sqrt 2}.
#' @inheritParams deconvolve
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
#' Camirand, F., Carroll, R.J., and Delaigle, A. (2018). Estimating the  
#' distribution of episodically consumed food measured with errors.  
#' \emph{Manuscript.} 
#'
#' Delaigle, A. and Gijbels, I. (2007). Frequent problems in calculating
#' integrals and optimizing objective functions: a case study in density
#' deconvolution. \emph{Statistics and Computing}, 17, 349 - 355.
#' 
#' Delaigle, A. and Hall, P. (2008). Using SIMEX for smoothing-parameter choice
#' in errors-in-variables problems. \emph{Journal of the American Statistical
#' Association}, 103, 481, 280-287
#' 
#' Delaigle, A. and Hall, P. (2016). Methodology for non-parametric 
#' deconvolution when the error distribution is unknown. \emph{Journal of the 
#' Royal Statistical Society: Series B (Statistical Methodology)}, 78, 1, 
#' 231-252.
#' 
#' Delaigle, A., Hall, P., and Meister, A. (2008). On Deconvolution with  
#' repeated measurements. \emph{Annals of Statistics}, 36, 665-685 
#' 
#' Fan,  J.,  and Truong,  Y. K. (1993),  Nonparametric Regression With Errors
#' in Variables,  \emph{The Annals of Statistics}.  21,  1900-1925.
#' 
#' @author Aurore Delaigle, Timothy Hyndman, Tianying Wang
#' 
#' @example man/examples/reg_deconvolve_eg.R
#'
#' @export

reg_deconvolve <- function(Y, 
                           W1, 
                           W2 = NULL, 
                           xx = seq(min(W1), max(W1), length.out = 100), 
                           errortype = NULL, 
                           sd_U = NULL, 
                           phiU = NULL, 
                           bw = NULL, 
                           rho = NULL, 
                           n_cores = NULL, 
                           kernel_type = c("default", "normal", "sinc"), 
                           seed = NULL,
                           use_alt_SIMEX_rep_opt = FALSE) {

    # Partial matching ---------------------------------------------------------
    dist_types <- c("normal", "laplace")
    if (!is.null(errortype)) {
        errortype <- dist_types[pmatch(tolower(errortype), dist_types)]
        if (is.na(errortype)) {
            stop("Please provide a valid errortype.")
        }
    }

    kernel_type <- match.arg(kernel_type)

    # Determine error type provided --------------------------------------------
    if (!is.null(phiU)) {
        errors <- "hom"
    } else if (!is.null(errortype)) {
        errors <- "hom"
    } else if (!is.null(W2)) {
        errors <- "rep"
    } else {
        errors<- "not_defined_properly"
    }

    # Augment W1 with W2 if provided along with phiU or sd_U and errortype
    if (errors == "hom" & !is.null(W2)) {
        W1 <- c(W1, W2)
        W2 <- NULL
        if (!is.null(phiU)) {
            warning("Both phiU and W2 have been provided. Continuing using errors defined by phiU and augmenting W1 with the data in W2.")
        } else {
            warning("Errortype and sd_U as well as W2 have been provided. Continuing using errors defined by errortype and sd_U and augmenting W1 with the data in W2.")
        }
    }

    if (!is.null(phiU) & !is.null(errortype) & !is.null(bw) & !is.null(rho)) {
        warning("Both phiU and errortype provided. Continuing ignoring errortype.")
    }

    # Check inputs -------------------------------------------------------------
    if (errors == 'not_defined_properly') {
        if (is.null(errortype) & is.null(phiU) & is.null(W2)) {
            stop("You must provide either errortype, phiU, or W2.")
        }
    }
    
    
    if (errors == "hom") {
        if (is.null(phiU) & !is.null(errortype) & is.null(sd_U)) {
            stop("You must provide sd_U along with errortype.")
        }

        if ((is.null(bw) | is.null(rho)) & (is.null(sd_U) | is.null(errortype))) {
            stop("If the bw and rho are not provided then you must provide errortype and sd_U.")
        }
    }

    if (errors == "rep"){
        if (!(length(W1) == length(W2))) {
            stop("W1 and W2 must be the same length.")
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
    RK <- kernel_list$rk
    muK2 <- kernel_list$muk2
    deltat <- tt[2] - tt[1]

    # Convert errortype to phiU ------------------------------------------------
    if (is.null(phiU) & errors == "hom"){
        phiU <- create_phiU(errors = "hom", errortype, sd_U)
    }

    # Calculate bandwidth ------------------------------------------------------
    if (is.null(bw) | is.null(rho)) {
        outcome_tmp <- bandwidth(W1 = W1, 
                                 W2 = W2,
                                 errortype = errortype,
                                 sd_U = sd_U,
                                 phiU = phiU, 
                                 Y = Y, 
                                 algorithm = "SIMEX", 
                                 n_cores = n_cores, 
                                 kernel_type = kernel_type, 
                                 seed = seed,
                                 use_alt_SIMEX_rep_opt = use_alt_SIMEX_rep_opt)
        bw <- outcome_tmp$h
        rho <- outcome_tmp$rho
    }

    # Compute estimate for m(X) ------------------------------------------------
    if (errors == "rep") {
        diff <- W1 - W2
        sd_U <- sqrt(stats::var(diff)/2)
        n <- length(W1)
        hnaive <- ((8 * sqrt(pi) * RK/3/muK2^2)^0.2) * 
            sqrt(stats::var(W1)) * n^(-1/5)
        h_min <- hnaive / 3
        t_search <- tt/h_min
        phi_U <- create_replicates_phi_U(W1, W2, t_search)

        if (use_alt_SIMEX_rep_opt) {
            W_bar <- (W1 + W2)/2
            y <- NWDecUknown(xx, W_bar, Y, phi_U, bw, rho, phiK, tt, deltat)
        } else {
            W_full <- c(W1, W2)
            Y_full <- c(Y, Y)
            y <- NWDecUknown(xx, W_full, Y_full, phi_U, bw, rho, phiK, tt, deltat)
        }
    } else {
        y <- NWDecUknown(xx, W1, Y, phiU, bw, rho, phiK, tt, deltat)    
    }

    structure(list(pdf = y, bw = bw, rho = rho, x = xx), 
                   class = c("reg_deconvolve", "list"))


}
