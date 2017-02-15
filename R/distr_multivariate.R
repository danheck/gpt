
# product-distributions for latent densities 

# Y: matrix with dimensions N x n.cont
# latent: list with univariate distributions (S4 class "d.uni")
# eta: unnamed (!) vector of parameter values
# const: named vector of parameter values
# @importFrom mvtnorm dmvnorm  ##  vcov=NULL, 

d.multi <- function(y, distr, eta, const, log = TRUE){
  
  lik <- mapply(dens, distr=distr, 
                y=data.frame(y),  #as.list(data.frame(as.matrix(y)))
                MoreArgs = list(eta=eta, const=const, log=log))
  if(nrow(y) == 1)
    lik <- matrix(lik, 1)
  
  # product distributions: multiply likelihood!
  if(log){
    ll <- rowSums(lik)
  }else{
    ll <- apply(lik, 1, prod)
  }
  return(ll)
}

# generate multivariate data:
r.multi <- function(n, distr, eta, const){
  Y <- mapply(rand, distr=distr, 
              MoreArgs = list(n=as.integer(n), eta=eta, const=const))
  return(Y)
}
