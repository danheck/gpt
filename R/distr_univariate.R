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



############## shifted LOGNORMAL ############ 
# mean, sd, shift
dlognorm <- function(x, mu = 4, sigma = 1, shift = 0, log=FALSE){
  if (log)
    suppressWarnings(dd <- dnorm(log(x - shift), mean = mu, sd = sigma, log=TRUE)- log(x - shift) )
  else
    suppressWarnings(dd <- dnorm(log(x - shift), mean = mu, sd = sigma) /(x - shift) )
  dd[is.na(dd)] <- ifelse(log, -Inf, 0)
  dd 
}
rlognorm <- function(n, mu = 4, sigma = 1, shift = 0){
  shift + rnorm(n, mean = mu, sd = sigma)
}
plognorm <- function(q, mu = 4, sigma = 1, shift = 0, log.p = FALSE){
  suppressWarnings(p <- pnorm(log(q - shift), mean = mu, sd = sigma, log.p = log.p))
  p[is.na(p)] <- ifelse(log.p, -Inf, 0)
  p
}
# curve(gpt:::dlognorm(x, 4.5, 1, 00), 0, 2000)
# curve(gpt:::plognorm(x, 4.5, 1, 00), 0, 2000)


########### shifted WEIBULL  ############ 
dsweibull <- function(x, shape, scale, shift, log = FALSE){
  dweibull(x-shift, shape=shape, scale=scale, log = log)
}
rsweibull <- function(n, shape, scale, shift){    
  shift + rweibull(n, shape=shape, scale=scale)
}
psweibull <- function(q, shape, scale, shift, log.p = FALSE){    
  pweibull(q - shift, shape=shape, scale=scale, log.p = log.p)
}
# curve(gpt:::psweibull(x, 1.5, 500, 300), 0, 2000)
# curve(gpt:::dsweibull(x, 1.5, 500, 300), 0, 2000)


######### shifted GAMMA
dsgamma <- function(x, shape, scale, shift, log = FALSE){
  dgamma(x - shift, shape=shape, scale=scale, log=log)
}
rsgamma <- function(n, shape, scale, shift){    
  shift + rgamma(n, shape = shape, scale = scale)
}
psgamma <- function(q, shape, scale, shift, log.p = FALSE){
  pgamma(q - shift, shape = shape, scale = scale, log = log.p)
}
# curve(gpt:::dsweibull(x, 1.5, 500, 300), 0, 2000)



