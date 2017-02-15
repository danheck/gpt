

############## LOGNORMAL ############ 
# mean, sd, shift
dlogNormal <- function(x, mu = 1, sigma=.1, nu=0, log=FALSE){
  suppressWarnings(
    dd <- dnorm(log(x-nu), 
                mean = mu, 
                sd = sigma, log=log)
  )
  dd[is.na(dd)] <- ifelse(log, -Inf, 0)
  dd
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


