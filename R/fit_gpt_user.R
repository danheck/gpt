####################################### Estimate MPT-continuous-mixture ########################

#' Fit GPT Model
#' 
#' @param x vector of observed discrete observations (e.g., choices) or name of variable in \code{data}
#' @param y vector  or matrix of observed continuous variable(s) (e.g., response times, confidence,...) or name of variable in \code{data}
#' @param data data frame containing the variables as named in \code{x} and \code{y}
#' @param file a character vector specifying the path to the model file
#' @param group vector with the number of observations as in \code{x} and \code{y}, indicating to which group each observation belongs
#' @param latent type of latent continuous distribution (one of \code{"normal"}, \code{"exgauss"}, \code{"exwald"}, \code{"gamma"}, \code{"weibull"}, \code{"lognormal"}, \code{"wald"}, or \code{beta}). Can be a  vector if multiple continuous variables \code{y} have different distributions (e.g., \code{latent = c("normal", "gamma")}) 
#' @param  restrictions list with parameter restrictions (e.g., \code{list("g=0.5", "mean1=mean2=100")})
#' @param baseline whether to fit baseline models that assume a saturated MPT structure and (a) separate continuous distributions per category (\code{y.per.cat} = no mixture) or (b) a single distribution across all categories (\code{y.null}). In each case, separate distributions are assumed for each continuous variable.
#' @param starting.values starting values for theta and eta (used only in first EM run). Note that sensible starting values are guessed by default.
#' @param eta.lower lower bound vector or scalar for eta parameters (assigned by name; or in alphabetical order). It is sufficient to constrain a subset of eta parameters selectively using a named vector.
#' @param eta.upper upper bound vector or scalar for eta parameters
#' @param EM.tol tolerance for EM algorithm
#' @param n.fit number of EM and optim fitting runs (if any \code{n.fit>1}: with random starting values)
#' @param maxit maximum number of EM and optim iterations, respectively
#' @param cpu number of cores used to fit individual data (default: number of cores minus one). Alternatively: a cluster spawned by \code{makeCluster()}
#' @param print whether to print current status
#' 
#' @details The following paramaters are used to specify the latent distributions (will be assigned in this order to the parameters listed in the model file as last argument in each row/branch):
#' \itemize{
#' \item{\code{"normal"}:}{mean (mu) and SD (sigma)}
#' \item{\code{"exgauss"} (sum of normal and independent exponential):}{mean and SD of normal (mu, sigma); mean of exponential (tau)}
#' \item{\code{"gamma"}:}{ shape (k), scale (theta), shift}
#' \item{\code{"weibull"}:}{shape (k), scale (lambda), shift}
#' \item{\code{"lognormal"}:}{mean and SD (of normal distribution before taking log), shift}
#' \item{\code{"wald"} (= inverse normal):}{mean (mu), shape (lambda), shift}
#' \item{\code{"exwald"} (= sum of inverse normal and independent exponential):}{mean, shape, shift}
#' \item{\code{"beta"}:}{shape1 (alpha) and shape2 (beta), see \link{dbeta}}
#' \item{\code{"vonmises"}:}{mu (mean), kappa (concentration) of von Mises distribution for circular data. }
#' }
#' @examples
#' \dontrun{
#' n <- c(targets=75, lures=75)     # number of items
#' theta <- c(do=.6,dn=.4, g=.5)          # MPT parameters
#' eta <- c(mu=400, sig=50, lambda_do=300, 
#'          lambda_go=500, lambda_gn=500, 
#'          lambda_dn=300)          # exGaussian parameters
#' 
#' file <- paste0(path.package("gpt"), "/models/2htm_exgauss.txt")
#' gen <- gen.gpt(n=n, theta=theta, eta=eta, latent="exgauss", file=file)
#' fit <- fit.gpt(x=gen$x, y=gen$y, latent="exgauss", file=file, 
#'                restrictions=list("do=dn", "lambda_do=lambda_dn", 
#'                                  "lambda_go=lambda_gn"))
#' fit
#' }
#' @export
fit.gpt <- function(x, y, data, file, latent, 
                    group=NULL, restrictions=NULL, baseline=FALSE,
                    starting.values=NULL, eta.lower=NULL, eta.upper=NULL, 
                    n.fit=c(3,1), maxit=c(300, 2000), EM.tol=.01, cpu=NULL, print=FALSE){
  
  # build model structure:
  gpt <- new("gpt", file, latent, restrictions)
  
  # check data input:
  data.checked <- data.check(mpt=gpt@mpt, x=x, y=y, 
                             data=data, group=group)
  x <- data.checked$x 
  y <- data.checked$y
  group <- data.checked$group
  
  t0 <- Sys.time()
  if(print) cat("Start fitting model: ", format(t0), "\n")
  # fit single data set:
  if(is.null(group) || length(unique(group))==1){
    # adjust bounds for eta parameters:
    gpt <- adjust.bounds(gpt = gpt, y = y, 
                         eta.lower = eta.lower, eta.upper = eta.upper) 
    
    fit <- fit.gpt.core(x=x,y=y,  gpt=gpt, 
                        baseline=baseline, starting.values=starting.values, 
                        n.fit=n.fit, EM.tol=EM.tol, maxit=maxit, print = print)
    
  }else{
    # multiple data sets:  
    fit.i <- function(i){
      sel <- group == unique(group)[i]
      xx <- x[sel]
      yy <- y[sel,,drop=FALSE]
      gpt <- adjust.bounds(gpt = gpt, y = yy, 
                           eta.lower = eta.lower, eta.upper = eta.upper) 
      
      fit.gpt.core(x=xx,y=yy, gpt = gpt, 
                   baseline=baseline, starting.values=starting.values, 
                   n.fit=n.fit, EM.tol=EM.tol, maxit=maxit)
    }
    
    if(missing(cpu) || is.null(cpu)){
      cpu <-  detectCores()-1
    }
    
    
    if(! "cluster" %in% class(cpu)){
      # number of CPUs provided as input:
      cl <- makeCluster(spec = cpu)
      stop.cl <- TRUE
    }else{
      # full cluster provided by input:
      cl <- cpu
      cpu <- length(cl)
      stop.cl <- FALSE
    }
    if(cpu>1){
      clusterExport(cl, varlist=c("x","y","group", "restrictions",
                                  "latent", "file", "starting.values",
                                  "eta.lower", "eta.upper", "n.fit", "EM.tol",
                                  "maxit"),
                    envir =environment())
      tmp <- clusterEvalQ(cl, library(gpt))
      
      fit.list <- parSapplyLB(cl, 1:length(unique(group)), fit.i, 
                              simplify = FALSE, USE.NAMES = FALSE)
      if(stop.cl) 
        stopCluster(cl = cl)
    }else{
      fit.list <- sapply(1:length(unique(group)), fit.i, 
                         simplify = FALSE, USE.NAMES = FALSE)
    }
    names(fit.list) <- unique(group)
    ests <- combine.gpt(fit.list, par.names = c(gpt@theta, gpt@eta))
    # dimnames(ests)[3] <- unique(group)
    
    ## aggregate parameters
    fit <- list(estimates = ests,
                individual.fits = fit.list)
  }
  t1 <- Sys.time()
  if(print) 
    cat("\nFinished: ", format(t1), 
        "\nRequired time: ", format(t1-t0))
  
  # store data and other model information:
  res <- c(fit, 
           list(data = list(x=x, y=y, group=group),
                gpt=gpt,
                input = list(file=file, latent=latent, restrictions=restrictions)))
  
  class(res) <- "gpt.fit"
  return(res)
}


combine.gpt <- function(fit.list, par.names){
  
  # subjects, parameters:
  S <- length(fit.list)
  P <- length(par.names)
  # P <- length(fit.list[[1]]$fit.EM$par)
  # par.names <- names(fit.list[[1]]$fit.EM$par)
  
  par <- array(NA, c(P, 2, S), dimnames=list(
    Parameter = par.names, 
    Statistic = c("Estimate", "SE"),
    Dataset = names(fit.list)))
  
  extract <- function(x, par, P){
    est <- x$fit.grad[[par]]
    if(is.null(est) || all(is.na(est))){
      return(rep(NA, P))
    }else{
      return(est)
    }}
  
  par[,1,] <- sapply(fit.list, extract, par="par", P=P)
  par[,2,] <- sapply(fit.list, extract, par="SE", P=P)
  return(par)
}

