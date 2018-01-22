kernel <- function(type = "default", coarse = FALSE){
  phik <- switch(type,
                 "default" = phiK2,
                 "normal" = normalK,
                 "sinc" = sincK)
  
  muk2 <- switch(type,
                 "default" = 6,
                 "normal" = 1,
                 "sinc" = NA)  # Integral does not converge

  rk <- switch(type,
               "default" = 1024 / 3003 / pi,
               "normal" = 1 / (2 * sqrt(pi)),
               "sinc" = 2)

  if (coarse) {
    tt <- switch(type,
                 "default" = seq(-1, 1, 1e-03),
                 "normal" = seq(-10, 10, length.out = 2000),
                 "sinc" = seq(-1, 1, 1e-03))
  } else {
    tt <- switch(type,
                 "default" = seq(-1, 1, 2e-04),
                 "normal" = seq(-10, 10, length.out = 10000),
                 "sinc" = seq(-1,1, 2e-04))
  }

  list(phik = phik, muk2 = muk2, rk = rk, tt = tt)
}

#Fourier transform of commonly used kernel function in deconvolution. (a default setting for the input)
phiK2 <- function(t) {
  y <- (1 - t^2)^3
  y[abs(t) > 1] <- 0
  y
}

normalK <- function(t) {
  y <- exp(-t^2/2)
  y
}

sincK <- function(t) {
  t*0 + 1
}
