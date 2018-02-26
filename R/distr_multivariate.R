# multivariate product-distributions for latent densities 

# Y: matrix with dimensions N x n.contin
# latent: list with univariate distributions (S4 class "contin")
# eta: unnamed (!) vector of parameter values
# const: named vector of parameter values
# @importFrom mvtnorm dmvnorm  ##  vcov=NULL, 

dmultivar <- function(y, distr, eta, const, log = TRUE){
  # check: dens(distr[[1]], y = c(y), eta = eta, const = const)
  lik <- mapply(dens, distr = distr, y = data.frame(y),
                MoreArgs = list(eta = eta, const = const, log=TRUE), 
                SIMPLIFY = FALSE)
  
  # independent/product distributions: multiply density / sum of log-density
  ll <- rowSums(do.call("cbind", lik))
  if (!log) ll <- exp(ll)
  ll
}

# random generation
rmultivar <- function(n, distr, eta, const){
  y <- mapply(rand, distr=distr, 
              MoreArgs = list(n=as.integer(n), eta=eta, const=const),
              SIMPLIFY = FALSE)
  do.call("cbind", y)
}

pmultivar <- function(y, distr, eta, const, log.p = TRUE){
  lik <- mapply(cdf, distr = distr, y = data.frame(y),
                MoreArgs = list(eta=eta, const=const, log.p = TRUE), 
                SIMPLIFY = FALSE)
  ll <- rowSums(do.call("cbind", lik))
  # independent/product distributions: multiply cdf / sum of log-cdf
  if (log.p) ll <- exp(ll)
  ll
}