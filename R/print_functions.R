
#' Print GPT Model
#' 
#' Either estimated parameters (single data set) or summary of estimates (if model was fitted for a sample of participants using \code{group})
#' 
#' @param x fitted model of \code{\link{fit.gpt}}
#' @param ci probability for confidence interval 
#' @param digits number of digits to be printed
#' @param ... ignored
#' 
#' @method print gpt.fit   
#' @export
print.gpt.fit <- function(x, ci=.95, digits=3,...){
  # single fit
  if(is.null(x$data$group)){
    zz <- qnorm(1- (1- ci)/2)
    tab <- cbind(
      # Estimate.EM=x$fit.EM$par, 
      Estimate=x$fit.grad$par, 
      # SE=x$fit.EM$SE, 
      SE=x$fit.grad$SE,
      CI.lower=x$fit.grad$par - zz*x$fit.grad$SE,
      CI.upper=x$fit.grad$par + zz*x$fit.grad$SE)
    print(round(tab, digits))
  }else{
    # individual estimates
    cat("Summary of individual estimates (of optim):\n")
    summ <- apply(x$estimates[,"Estimate",], 1, summary, na.rm=TRUE)
    print(round(summ, digits = digits))
  }
}


# print.gpt <- function(x, ...){
#   print(x)
# }

# @importFrom base print
#  @export
setMethod(
  f="print",
  signature = "gpt",
  definition = function(x, ...){
    cat('Free MPT Parameters:\n')
    print(x@theta)
    cat('\nFree continuous Parameters:\n')
    print(x@eta)
    cat("\nLatent continuous structure: [", 
        paste(sapply(x@distr[[1]], 
                     function(dd) dd@label), collapse=" / "), 
        "], \n  with following parameterization:\n\n")
    X <- get.X(x@distr, x@eta, x@const)
    print(X)
  }
)


# get clean assignment structure for printing:
get.X <- function(distr, eta.names, const){
  
  # recover eta indices:
  X <- X.tmp <- t(sapply(distr, function(dd) 
    unlist(sapply(dd, function(d) d@eta.idx))))
  
  # replace by labels/constants:
  X[X.tmp>0] <- eta.names[ X.tmp[X.tmp>0]]
  X[X.tmp<0] <- const[-X.tmp[X.tmp<0]]
  
  # meaningful names:
  par.per.lat <- sapply(distr[[1]], 
                        function(dd) length(dd@eta.idx))
  colnames(X) <-  rep(paste0("c", 1:length(distr[[1]]), "_",
                             sapply(distr[[1]], function(dd) dd@label)) , 
                      par.per.lat)
  X
}