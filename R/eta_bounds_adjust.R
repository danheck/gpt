

# adjust eta boundaries after the data are known and checked
adjust.bounds <- function(gpt, y, eta.lower = NULL, eta.upper = NULL){
  
  npar <- length(gpt@eta)
  
  if(!is.null(eta.lower) && length(eta.lower) == 1 && is.null(names(eta.lower))){
    eta.lower <- rep(eta.lower, npar)
  }
  if(!is.null(eta.upper) && length(eta.upper) == 1 && is.null(names(eta.upper))){
    eta.upper <- rep(eta.upper, npar)
  }
  
  # check parameter names etc.
  eta.lower <- check.input.par(par = eta.lower, names = gpt@eta, replace.NA = -Inf)
  eta.upper <- check.input.par(par = eta.upper, names = gpt@eta, replace.NA = Inf)
  
  orderOK <- all(eta.upper > eta.lower, na.rm = TRUE)
  if(!is.null(eta.lower) && !is.null(eta.upper) && !orderOK)
    stop("some values of eta.lower are larger or equal to those of eta.upper!")
  
  
  
  for (ss in seq_along(gpt@distr)){
    for (cc in seq_along(gpt@distr[[ss]])){
      tmp <- gpt@distr[[ss]][[cc]]
      
      
      # boundaries due to univariate distribution parameters:
      gpt@distr[[ss]][[cc]] <-  make.distr(label = tmp@label, 
                                           eta.idx = tmp@eta.idx,
                                           y = y[,cc])
      
      pars <- gpt@eta.repar[tmp@eta.idx]
      
      # check whether user-specified bounds are more strict than necessary default boundaries:
      if (!is.null(eta.lower) && length(eta.lower) > 0 && !anyNA(eta.lower[pars])){
        gpt@distr[[ss]][[cc]]@lower <- 
          pmax(eta.lower[pars], 
               gpt@distr[[ss]][[cc]]@lower, na.rm = TRUE)
      }
      
      if (!is.null(eta.upper)&& length(eta.upper) > 0 && !anyNA(eta.upper[pars])){
        gpt@distr[[ss]][[cc]]@upper <- 
          pmin(eta.upper[pars], 
               gpt@distr[[ss]][[cc]]@upper, na.rm = TRUE)
      }
    }
  }
  
  gpt
}


