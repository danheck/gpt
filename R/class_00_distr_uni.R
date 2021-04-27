setClass("contin", 
         representation(label = "character", 
                        dens = "character",
                        eta.idx = "integer",
                        range = "numeric",
                        lower = "numeric",
                        upper = "numeric"),
         
         prototype(label = NA_character_, 
                   dens = NA_character_, 
                   eta.idx = NA_integer_,
                   range = c(-Inf, Inf),
                   lower = NA_real_,
                   upper = NA_real_),
         
         validity=function(object){
           
           if(length(object@range) != 2)
             return("Range requires two values!")
           
           if(any(object@lower >= object@upper))
             return("Lower bounds for parameters are not smaller than upper bounds!")
           
           
           ee <- length(object@eta.idx)
           
           if(ee != length(object@lower) || ee != length(object@upper))
             return("Lengths of eta.idx, lower and upper bounds do not match!")
           
           if(any(is.na(object@lower), is.na(object@upper)))
             return("Lower and upper bounds missing or not valid numbers.")
           
           return(TRUE)
         })


################# DISTRIBUTIONS

#' @importFrom stats dnorm rnorm dbeta rbeta dgamma rgamma dweibull rweibull runif dunif rexp dexp integrate na.omit pchisq pnorm quantile
# define distribution object with correct density and range:
make.distr <- function(label, eta.idx, y=NULL){
  
  # for shifted distributions:
  if(is.null(y)){
    miny <- Inf
  }else{
    miny <- min(y)-.001
  }
  switch(label,
         "normal" = {
           dens <- "norm"
           range <- c(-Inf, Inf)
           lower <- c(mean = -Inf, sd = 0)
           upper <- c(mean = Inf, sd = Inf)},
         "lognormal" = {
           dens <- "lognorm" 
           range <- c(0, Inf)
           lower <- c(mu = -Inf, sigma = 0, shift = 0)
           upper <- c(mu = Inf, sigma = Inf, shift = miny)},    # shift range!
         "beta" = {
           dens <- "beta"  
           range <- c(0, 1)
           lower <- c(shape1 = 0, shape2 = 0)
           upper <- c(shape1 = Inf, shape2 = Inf)},
         "gamma" = {
           dens <- "sgamma"
           range <- c(0, Inf)
           lower <- c(shape = 0, scale = 1e-8, shift = 0)     # shift range!
           upper <- c(shape = Inf, scale = Inf, shift = miny)},
         "exgauss" = {
           dens <- "exgauss"
           range <- c(0, Inf)
           lower <- c(mu = -Inf, sigma = 1e-8, tau = 1e-3)
           upper <- c(mu = Inf, sigma = Inf, tau = Inf)},
         "wald" = {
           dens <- "wald"
           range <- c(-Inf, Inf)
           lower <- c(mean = 1e-8, shape = 1e-8, shift = 0)
           upper <- c(mean = Inf,  shape = Inf, shift = miny)},         # shift range!
         "exwald" = {
           dens <- "exwald"
           range <- c(0, Inf)
           lower <- c(m = 1e-8, a = 1e-8, t = 1e-3)
           upper <- c(m = Inf, a = Inf, t = Inf)},
         "weibull" = {
           dens <- "sweibull"
           range <- c(0, Inf)
           lower <- c(shape = 1e-8, scale = 1e-8, shift = 0)
           upper <- c(shape = Inf, scale = Inf, shift = miny)},       # shift range!
         "mises" = {
           dens <- "mises"
           range <- c(-pi, pi)
           lower <- c(mu = -pi, kappa = 0)
           upper <- c(mu = pi, kappa = Inf)},       # shift range!
         "unif" = {
           dens <- "unif"
           range <- c(-Inf, Inf)
           lower <- c(min = -Inf, max = -Inf)
           upper <- c(min = Inf, max = Inf)},
         {dens <-  label
         range <- c(-Inf, Inf)
         lower <- rep(-Inf, length(eta.idx))
         upper <- rep(Inf, length(eta.idx))
         # warning ("Distribution label not recognized.", 
         #          "\n Trying to use  'd", dens, "'  as density and",
         #          "\n 'd", dens, "' for data generation with three parameters.")
         }
  )
  new("contin", label=label, dens=dens, eta.idx=eta.idx, 
      range=range, lower=lower, upper=upper)
}



# setMethod(
#   f="identical",
#   signature="distr.uni",
#   definition = function(x, y){
#     ident <- 
#       x@label == y@label &&
#       x@dens == y@dens &&
#       x@eta.idx == y@eta.idx &&
#       x@range == y@range &&
#       x@lower == y@lower &&
#       x@upper == y@upper
#     return(ident)
#   }
# )
