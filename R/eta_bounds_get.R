
# NOTE: user.defined overwrites the automatic search for boundaries (necessary for reparameterization)
get.eta.bound <- function (gpt, lower = TRUE, user.defined = NULL, warning = FALSE){
  
  bnd <- rep(ifelse(lower, -Inf, Inf), length(gpt@eta))
  names(bnd) <- gpt@eta
  
  if (!all(gpt@eta.repar %in% gpt@eta)){
    # only show warning once:
    if (warning)
      warning("Please make sure to define sensible bounds via 'eta.lower'/'eta.upper' \n",
              "  when using reparameterized functions for the component parameters eta!")
    if (!missing(user.defined) && !is.null(user.defined))
      bnd[names(user.defined)] <- user.defined #[gpt@eta]
    
  } else {
    for(base in gpt@distr){
      for(cc in base){
        eta.free <- cc@eta.idx > 0
        if (any(eta.free)){
          # eta.idx <- cc@eta.idx[eta.free]  # old version without eta.repar
          eta.idx <- gpt@eta.repar[cc@eta.idx[eta.free]]
          if (lower)
            bnd[eta.idx] <- pmax(bnd[eta.idx], cc@lower[eta.free])
          else
            bnd[eta.idx] <- pmin(bnd[eta.idx], cc@upper[eta.free])
        }
      }
    }
  }
  bnd
}


# get vector of lower bounds for gpt model
get.eta.lower <- function(gpt){
  
  lower <- rep(-Inf, length(gpt@eta))
  for(base in gpt@distr){
    for(cc in base){
      eta.free <- cc@eta.idx > 0
      if (any(eta.free)){
        eta.idx <- cc@eta.idx[eta.free]
        lower[eta.idx] <- pmax(lower[eta.idx], cc@lower[eta.free])
      }
    }
  }
  # extract bounds:
  # mat <- sapply(gpt@distr,
  #               function(dlist) sapply(dlist, 
  #                                      function(dd)
  #                                        dd@lower))
  # # get indices:
  # idx <- sapply(gpt@distr,
  #               function(dlist) sapply(dlist, 
  #                                      function(dd)
  #                                        dd@eta.idx))
  # # use maximum if bounds differ:
  # lower <- as.vector(by(c(mat)[c(idx)>0], 
  #                       c(idx)[c(idx)>0], max))
  lower
}


get.eta.upper <- function(gpt){
  
  upper <- rep(Inf, length(gpt@eta))
  for(base in gpt@distr){
    for(cc in base){
      eta.free <- cc@eta.idx > 0
      if (any(eta.free)){
        eta.idx <- cc@eta.idx[eta.free]
        upper[eta.idx] <- pmin(upper[eta.idx], cc@upper[eta.free])
      }
      # upper[cc@eta.idx] <- pmax(upper[cc@eta.idx], cc@upper)
    }
  }
  # extract bounds:
  # mat <- sapply(gpt@distr,
  #               function(dlist) sapply(dlist, 
  #                                      function(dd)
  #                                        dd@upper))
  # # get indices:
  # idx <- sapply(gpt@distr,
  #               function(dlist) sapply(dlist, 
  #                                      function(dd)
  #                                        dd@eta.idx))
  # # use maximum if bounds differ:
  # upper <-  as.vector(by(c(mat)[c(idx)>0], 
  #                        c(idx)[c(idx)>0], min))
  upper
}