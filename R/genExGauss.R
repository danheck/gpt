# 
# ## special densities
# 
# # Generalized ex-Gaussian
# # 
# # Distribution of the sum of one normal and several exponential random variables (all independent).
# # @param x, q vector of quantiles
# # @param n number of observations
# # @param mean mean of normal component
# # @param sd standard deviation of normal component
# # @param lambda vecor of means of one or more exponential components
# # @param log if TRUE, log densities are returned
# # @param prec precision of approximation of density: minimal difference between two values of the vector \code{lambda}
# # @export
# dGenExGauss <- function(x, mean, sd, lambda, log=F, prec=1e-5){
#   
#   lambda <- lambda[lambda!=0]
#   
#   
#   # adjsut lambda: precision of approximation
#   L <- length(lambda)
#   differ <- dist(lambda) 
#   if (L>1 && min(differ)< prec ){
#     selmat <- as.matrix(differ) < prec
#     selmat[upper.tri(selmat, diag=T)] <- F
#     idx <- apply(selmat, 1, any)
#     lambda[idx] <- lambda[idx] + prec*(1:sum(idx))
#   }
#   
#   lambda <- 1/lambda
#  
#   d <- 0
#   for (i in 1:L){
#     e <- prod(1/(lambda[-i]-lambda[i]))
#     f <- pnorm(x, mean+lambda[i]*sd^2, sd, log.p = T)
#     g <- lambda[i]*(mean-x)+lambda[i]^2*sd^2/2
#     d <- d +  e*exp(g+f)
#   }
#   
#   if (log) {
#     sum(log(lambda)) +log(d)
#   }else{ 
#     prod(lambda)*d
#   }
# }
# 
# # @param q vector of quantiles
# # @param lower lower bound for numerical integration of density expressed as a distance in standard deviation from the mean of te normal component: mean-lower*sd
# # @rdname  dGenExGauss
# # @export
# pGenExGauss <- function(q, mean, sd, lambda, prec=1e-5,lower=5){
#   p <- sapply(q, FUN=function(qq) integrate(dGenExGauss, lower=mean-lower*sd,upper=qq, 
#                                            mean=mean, sd=sd, lambda=lambda)$value)
#   p[p<=0] <- 0
#   p[p>=1] <- 1
#   return(p)
# }
# 
# # @rdname  dGenExGauss
# # @export
# rGenExGauss <- function(n, mean, sd, lambda){
#   if(!is.matrix(lambda)){
#     lambda <- matrix(lambda, n, length(lambda), byrow=TRUE)
#   }
# 
#   exps <- apply(lambda, 2, function(ll) rexp(n, 1/ll))
#   if(n == 1){
#     x <- rnorm(n, mean, sd) + sum(exps)
#   }else{
#     x <- rnorm(n, mean, sd) + rowSums(exps)
#   }
#   return(x)
# }
