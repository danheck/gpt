
# get vector of lower bounds for gpt model
get.eta.lower <- function(gpt){
  
  lower <- rep(-Inf, length(gpt@eta))
  for(base in gpt@distr){
    for(cc in base){
      lower[cc@eta.idx] <- pmax(lower[cc@eta.idx], cc@lower)
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
      upper[cc@eta.idx] <- pmax(upper[cc@eta.idx], cc@upper)
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