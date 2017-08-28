############## von Mises ############ 
# mean, sd, shift
dmises <- function(x, mu = 0, kappa=1, log=FALSE){
  if (kappa < 0){
    stop("Parameter 'kappa' must be nonnegative.")
  } else if (kappa == 0){
    dd <- rep(-log(2*pi), length(x))
  } else {
    suppressWarnings(
      dd <- kappa * cos(x - mu) - log(2*pi) - log(besselI(kappa, mu))
    )
  }
  if (!log)
    dd <- exp(dd)
  dd[is.na(dd)] <- ifelse(log, -Inf, 0)
  return(dd)
}

# Best, D. J., & Fisher, N. I. (1979). Efficient Simulation of the von Mises Distribution. Journal of the Royal Statistical Society. Series C (Applied Statistics), 28(2), 152-157.
rmises <- function(n, mu = 0, kappa=1){
  if (kappa < 0){
    stop("Parameter 'kappa' must be nonnegative.")
  } else if (kappa == 0){
    x <- runif(n, -pi, pi)
    return(x)
  } else{
    x <- rep(NA, n)
    a <- 1 + (1 + 4 * (kappa^2))^0.5
    b <- (a - (2 * a)^0.5)/(2 * kappa)
    r <- (1 + b^2)/(2 * b)
    idx <- 1:n
    while ( length(idx) > 0) {
      U1 <- runif(length(idx))
      z <- cos(pi * U1)
      f <- (1 + r * z)/(r + z)
      c <- kappa * (r - f)
      U2 <- runif(length(idx))
      sel <- (c * (2 - c) - U2 > 0) | (log(c/U2) + 1 - c >= 0 )
      if (any(sel)){
        U3 <- runif(sum(sel))
        x[idx[sel]] <- sign(U3 - 0.5) * acos(f[sel]) + mu + pi
      }
      idx <- which(is.na(x))
    }
    x <- x %% (2 * pi)
    return (x - pi)
  }
}



############## LOGNORMAL ############ 
# mean, sd, shift
dlogNormal <- function(x, mu = 1, sigma=.1, nu=0, log=FALSE){
  suppressWarnings(
    dd <- dnorm(log(x-nu), 
                mean = mu, 
                sd = sigma, log=log)
  )
  dd[is.na(dd)] <- ifelse(log, -Inf, 0)
  return(dd)
}
rlogNormal <- function(n, mu = 1, sigma=.1, nu=0){
  x <- rnorm(n, mean=mu, sd=sigma)
  return(nu + exp(x))
}


########### WEIBULL  ############ 
dWeibull <- function(x, shape, scale, shift, log=FALSE){
  lik <- dweibull(x-shift, 
                  shape=shape, 
                  scale=scale, log=log)
  return(lik)
}
rWeibull <- function(n, shape, scale, shift){    
  x <- shift + rweibull(n, shape=shape, scale=scale)
  return(x)
}

######### shifted gamma
dsgamma <- function(x, shape, rate, shift, log=FALSE){
  lik <- dgamma(x-shift, 
                shape=shape,
                rate=rate, 
                log=log)
  return(lik)
}
rsgamma <- function(n, shape, rate, shift){    
  x <- shift + rgamma(n, shape=shape, rate=rate)
  return(x)
}


