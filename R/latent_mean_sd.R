#' Mean and SD of Basis Distributions
#' 
#' Computes summary statistics for basis distributions, conditional on the 
#' parameter estimates.
#' 
#' @param model a fitted GPT model. See \code{\link{fit.gpt}}
#' @param base optional: a vector with indices for component distributions
#' @param dim dimension of the continuous variable 
#' @export
base_mean_sd <- function(model, base, dim = 1){
  
  distr <- model$gpt@distr
  if (missing(base)) 
    base <- seq.int(length(distr))
  B <- length(base)
  stats <- matrix(NA, B, 2, 
                  dimnames = list("base" = seq.int(B), 
                                  "statistic" = c("mean", "sd"))) #, "SE(mean)", "SE(sd)")))
  eta <- model$fit.grad$par[model$gpt@eta]
  for (i in seq.int(B)){
    # get parameters in branch i
    stats[i,] <- latent.mean.sd(distr[[i]][[dim]], eta, model$gpt@const)
    
  }
  stats
}

# transform mean/SD + SE
latent.mean.sd <- function(distr.uni, eta, const){
  
  # get parameters of fitted model
  names(eta) <- names(const) <- NULL
  par <- distr.uni@eta.idx
  cc <- distr.uni@eta.idx < 0
  # free parameters:
  par[!cc] <- eta[distr.uni@eta.idx[!cc]]
  # constants:
  par[cc] <- const[-distr.uni@eta.idx[cc]]

  ms <- switch(distr.uni@label,
               "normal" = {
                 c("mean" = par[1], "sd" = par[2])
               },
               "lognormal" = {
                 mu <- par[1] ; sigma <- par[2] ; shift <- par[3]
                 c("mean" = exp(mu + sigma^2/2) + shift,
                   "SD" = sqrt( (exp(sigma^2) - 1) * exp(2*mu + sigma^2)))
               },
               "beta" = {  
                 s1 <- par[1] ; s2 <- par[2]  # wiki: alpha/beta
                 c("mean" = s1/(s1 + s2),
                   "sd" = sqrt(s1*s2 / ((s1+s2)^2 * (s1+s2+1))))
              },
               "gamma" = {  # wiki: k = shape, theta = scale
                 shape <- par[1] ; scale <- par[2] ; shift <- par[3] 
                 c("mean" = shape * scale + shift,
                   "sd" = sqrt(shape * scale^2))
               },
               "exgauss" = {
                 mu <- par[1] ; sigma <- par[2] ; tau <- par[3]
                 c("mean" = mu + tau,
                   "sd" = sqrt(sigma^2 + tau^2))
               },
               "wald" = {
                 mean <- par[1] ; shape <- par[2] ; shift <- par[3]
                 c("mean" = mean + shift,
                   "sd" = sqrt(mean^3 / shape))
               },
               "exwald" = {
                 m <- par[1]
                 a <- par[2]
                 t <- par[3]
                 c("mean" = t + a/m,
                   "sd" = sqrt(t^2 + a/m^3))
               },
               "weibull" = {  # wiki: k=shape , lambda = scale
                 shape <- par[1] ; scale <- par[2] ; shift <- par[3]
                 c("mean" = scale * gamma(1 + 1/scale) + shift,
                   "sd" = sqrt(scale^2 * (gamma(1+2/shape) - gamma(1+1/shape)^2)))
               },       # shift range!
               "mises" = {
                 mu <- par[1] ; kappa <- par[2]
                 c("mean" = mu, 
                   "sd" = sqrt(1 - besselI(kappa, 1)/besselI(kappa, 0)))
               },
               "unif" = {
                 min <- par[1] ; max <- par[2]
                 c("mean" = (max + min)/2,
                   "sd" = sqrt((max - min)^2 / 12))
               },
               {mean <- sd <- NA})
  
  names(ms) <- c("mean", "sd") #, "SE.mean", "SE.sd")
  ms
}