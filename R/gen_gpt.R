
#' Generate Single Data Set
#' 
#' @inheritParams fit.gpt
#' @param n vector of the length of number of trees with n per tree
#' @param theta vector of true MPT parameter values (consider also restricted parameters!). Values will be appropriatly reordered if a named vector ist given (otherwise, check order!).
#' @param eta vector of true continuous parameters. Values will be appropriatly reordered if a named vector ist given (otherwise, check order!).
#' @examples
#' ###### 2-High-Threshold Model (+2 normal distributions) ######
#' 
#' ### parameters
#' n <- c(targets=20, lures=20)     # number of items
#' theta <- c(do=.6, g=.5)          # MPT parameters
#' eta <- c(m1_d=100, m1_g=110, s1=10,
#'          m2_d=30, m2_g=50, s2=5)     # normal distributions
#' file <- paste0(path.package("gpt"), "/models/2htm_2normal.txt")
#' 
#' gen <- gen.gpt(n=n, theta=theta, eta=eta, file=file,
#'               latent=c("normal", "normal"), restrictions=list("do=dn"))
#' head(gen) 
#' # check means of latent continuous distributions:
#' tapply(gen$y.1, gen$state, mean)
#' 
#' @seealso gen.gpt.sample
#' @export
gen.gpt <- function(n, theta, eta, file, latent, restrictions=NULL){
  
  # build S4 model:
  gpt <- new("gpt", file=file, latent=latent, 
                restrictions=restrictions)
  
  # check parameters
  theta.names <- names(gpt@mpt@theta)[gpt@mpt@theta == -.5]
  theta <- check.input.par(theta, theta.names)
  eta <- check.input.par(eta, gpt@eta)
  
  # use S4 data generation for gpt:
  data <- rand(distr=gpt, n=n, eta=eta, theta=theta)
  return(data)
}


