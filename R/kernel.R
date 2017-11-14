kernel <- function(type = "default", coarse = FALSE){
  phik <- switch(type,
                 "default" = phiK2)
  
  muk2 <- switch(type,
                 "default" = 6)

  rk <- switch(type,
               "default" = 1024 / 3003 / pi)

  if (coarse) {
    tt <- switch(type,
                 "default" = seq(-1, 1, 1e-03))
  } else {
    tt <- switch(type,
                 "default" = seq(-1, 1, 2e-04))
  }

  list(phik = phik, muk2 = muk2, rk = rk, tt = tt)
}

#Fourier transform of commonly used kernel function in deconvolution. (a default setting for the input)
phiK2 <- function(t) {
  y <- (1 - t^2)^3
  y[abs(t) > 1] <- 0
  y
}
