

##### 3.  continuous structure

make.gpt <- function(tab,  latent){
  
  distr <- c("exgauss", "exwald", "gamma", "weibull", 
             "lognormal", "wald", "normal", "beta", "mises", "unif")
  if(any(! latent %in% distr))
    warning(paste("Check definition of latent distribution functions!",
                  "\n  Should be one of those listed in help(fit.gpt).",
                  "\n  Currently defined:",paste(latent,collapse=", ")))
  
  # latent: vector of latent distributions ("mvnormal" are collapsed including covariances)
  
  
  n.cont <- ncol(tab) - 3
  if(length(latent) == 1 ){
    latent <- rep(latent, n.cont)
  }else if (length(latent) != n.cont ){
    stop("Number of latent distributions does not match with the model file!")
  }
  
  B <- nrow(tab)
  
  # assign parameters to latent states:
  par.per.lat <- ifelse(latent %in% c("normal", "beta", "mises", "unif"), 2, 3)
  sel.cont <- rep(1:length(latent), par.per.lat)
  X.named <- X.full <- matrix(NA, B, ncol=sum(par.per.lat))
  for(cc in 1:n.cont){
    parlist <- strsplit(tab[,3+cc], ",")
    for(pp in 1:B){
      if(length(parlist[[pp]]) != par.per.lat[cc])
        stop("Number of parameters in model file (=",length(parlist[[pp]]), ")\n",
             "  does not match number of parameters of 'latent' distribution (=",par.per.lat[cc],").")
      X.named[pp,sel.cont == cc] <- parlist[[pp]]
    }
  }
  
  eta.names <- sort(na.omit(unique(c(X.named))))
  X.full <- matrix(match(X.named, eta.names), B, 
                   sum(par.per.lat))
  X.full[is.na(X.full)] <- 0
  # rownames(X.full) <- rownames(X.named) <- branch.names
  colnames(X.full) <- colnames(X.named) <- sel.cont
  
  unique.state <- !duplicated(X.full)
  if(any(X.full== 0) ){
    stop(paste("  Check number of latent continuous parameters in model file!",
               "\n  Should be [", paste(par.per.lat,collapse=" / "),"] per line for [",
               paste(latent, collapse=" / "), "]"))
  }
  
  # reduce duplicated states (for speed)
  X <- X.full[unique.state,,drop=FALSE]   # can be recovered by get.X()
  
  # number of distinct basis distributions:
  S <- nrow(X)
  map.vec <- rep(NA, B)
  for(bb in 1:B){
    sel.state <- apply(X, 1, function(xx) all(xx == X.full[bb,]))
    map.vec[bb] <- (1:S)[sel.state]
  }
  
  # list of S4 distributions:
  distr <- vector("list", S)
  names(distr) <- paste0("base", 1:S)
  for(s in 1:S){
    distr[[s]] <- vector("list", n.cont)
    names(distr[[s]]) <- paste0("cont", 1:n.cont)
    for(cc in 1:n.cont){
      distr[[s]][[cc]] <- make.distr(label = latent[cc], 
                                     eta.idx = as.integer(X[s, colnames(X) == cc]))
    }
  }
  
  # SxB matrix:
  # B-dim vector of states z_jn -> S-dim vector of latent continuous distributions
  # map <- matrix(0, S, B, dimnames=list(paste0("base", 1:S), paste0("br",1:B)))
  # for(bb in 1:B)
  #   map[map.vec[bb], bb] <- 1
  
  res <- list(
    # map branches to basis distributions:
    map.vec=map.vec, 
    # distributions:
    distr=distr, 
    # labels:
    eta.names=eta.names)
  res
}

