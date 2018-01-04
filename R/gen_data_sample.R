
#' Generate Sample of Data Sets
#' 
#' Allows to simulate parameter heterogeneity. Note that truncated normal distributions are used to ensure that parameters within the specified bounds.
#' 
#' @param S sample size (number of participants)
#' @param theta_sd vector giving the standard deviation of normally distributed MPT parameters across participants (default: homogeneity)
#' @param eta_sd vector giving the standard deviation of normally distributed continuous parameters across participants (default: homogeneity)
# @param eta.lower lower bounds for continuos parameters (might lead to a decrease in variance!)
# @param eta.upper upper bounds for continuous parameters
#' @param cpu number of cores used for data generation (default: number of cores minus one). alternatively, a parallel cluster initialized by `cl <- parallel::makeCluster(4)` 
#' @inheritParams gpt_fit
#' @inheritParams gpt_gen
#' 
#' @seealso gpt_gen
#' @examples
#' ###### 2-High-Threshold Model (with fixed guessing) ######
#' 
#' \dontrun{
#' n <- c(targets=20, lures=20)     # number of items
#' theta <- c(do=.6,dn=.4, g=.5)          # MPT parameters
#' eta <- c(mu=400, sig=50, lambda_do=200, 
#'          lambda_go=600, lambda_gn=400, 
#'          lambda_dn=300)          # exGaussian parameters
#' file <- paste0(path.package("gpt"), "/models/2htm_exgauss.txt")
#' 
#' gen <- gpt_gen_sample(S=3, n=n, theta=theta, eta=eta, 
#'                       theta_sd=.1, eta_sd=10,
#'                       file =file, latent="exgauss")
#' sapply(gen, head, 3) 
#' # check mean of latent continuous distributions:
#' by(gen$data$y, gen$data$state, mean)
#' }
#' @export
#' @import parallel
gpt_gen_sample <- function(S, n, theta, eta, theta_sd=0, eta_sd=0, 
                           file, latent, restrictions=NULL, cpu){
  
  # build S4 model:
  gpt <- new("gpt", file=file, latent=latent, 
                restrictions=restrictions)
  
  # check parameters
  theta.names <- names(gpt@mpt@theta)[gpt@mpt@theta == -.5]
  theta <- check.input.par(theta, theta.names)
  eta <- check.input.par(eta, gpt@eta)
  theta_sd <- check.input.par(theta_sd, theta.names)
  eta_sd <- check.input.par(eta_sd, gpt@eta)
  
  
  # generate true values with heterogeneity
  theta.sample <- gen.pars(S, theta, theta_sd, 
                           lower = rep(0,length(theta)), 
                           upper = rep(1, length(theta)))
  eta.sample <- gen.pars(S, eta, eta_sd, 
                         lower = get.eta.lower(gpt), 
                         upper = get.eta.upper(gpt))
  
  
  # function to generate data for single participant
  gen.i <- function(i){
    gen <- gpt_gen(n=n, theta=theta.sample[i,], eta=eta.sample[i,], 
                   file=file, latent=latent, restrictions=restrictions)
    data.frame(id = i, gen)
  }
  
  # multicore management
  if(missing(cpu) || is.null(cpu)) 
    cpu <- detectCores()-1
  
  if(! "cluster" %in% class(cpu)){
    cl <- makeCluster(spec = cpu)
    stop.cl <- TRUE
  }else{
    cl <- cpu
    cpu <- length(cl)
    stop.cl <- FALSE
  }
  if(cpu>1){
    clusterExport(cl, varlist=c("n","theta.sample","eta.sample", 
                                "latent", "file", "gen.i"),
                  envir =environment())
    tmp <- clusterEvalQ(cl, library(gpt))
    data.list <- parSapply(cl, 1:S, gen.i, simplify = FALSE, USE.NAMES = FALSE)
    if(stop.cl) stopCluster(cl = cl)
    
  }else{
    data.list <- sapply(1:S, gen.i, simplify = FALSE, USE.NAMES = FALSE)
  }
  
  data <- do.call("rbind", data.list)
  res <- list(data = data, theta = theta.sample, eta = eta.sample)
  return(res)
}



# generated set of normally distributed parameters
# (truncating based on necessary bounds for parameters)
gen.pars <- function(S, par, sd, lower, upper){
  
  if(any(lower >= upper)){
    warning("Lower bounds are equal or larger than upper bounds!")
  }
  sel.l <- par < lower 
  sel.u <- par > upper
  if(any(sel.l | sel.u)){
    warning("Mean for parameter distribution is outside of the boundaries!")
    if(any(sd[sel.l | sel.u] == 0)){
      warning("Unable to sample parameters from truncated normal distribution if",
              "\n  parameter mean is outside of the boundaries and SD=0.",
              "\n  The mean is shifted to the boundary automatically.")
      par[sel.l & sd[sel.l]==0] <- lower[sel.l & sd[sel.l]==0]
      par[sel.u & sd[sel.u]==0] <- upper[sel.u & sd[sel.u]==0]
    }
  }
  
  # initialize
  sample <- matrix(NA, S, length(par), 
                   dimnames=list(NULL, names(par)))
  
  # truncation: sample until values fulfill bounds
  while(any(is.na(sample))){
    
    # hierarchical distribution: normal
    tmp <- matrix(rnorm(S*length(par), par, sd), S, byrow=TRUE)
    
    # truncation:
    sel <- t(apply(tmp, 1, function(pp) pp >= lower & pp <= upper))
    sample[sel] <- tmp[sel]
  }
  
  sample
}

