

##### 3.  continuous structure

make.gpt <- function(tab, latent, restrictions = NULL){
  
  distributions <- c("exgauss", "exwald", "gamma", "weibull", 
                     "lognormal", "wald", "normal", "beta", "mises", "unif", "custom")
  if(any(! latent %in% distributions))
    warning(paste("Check definition of latent distribution functions!",
                  "\n  Should be one of those listed in help(gpt_fit).",
                  "\n  Currently defined:",paste(latent,collapse=", ")))
  
  # latent: vector of latent distributions ("mvnormal" are collapsed including covariances)
  
  
  
  #################### constraints
  # old coding:
  # 1,2,...    = free parameter:        eta[eta.idx]
  # -1,-2,...  = constant parameter:    const[-eta.idx]
  # new coding:
  # constants are replaced in the eta.repar
  if (!missing(restrictions) && !is.null(restrictions)){
    for (i in 1:length(restrictions)){
      restrictions[[i]] <- gsub(" ", "", restrictions[[i]], fixed = TRUE)
      vec <- unlist(strsplit(restrictions[[i]], "="))
      for (k in 1:(length(vec) - 1)){
        # replaced by last value of equalities: x=y=z
        tab$c1 <- gsub(paste0("\\b", vec[k], "\\b"), # regex: boundary of words (not points!!)
                       vec[length(vec)],
                       tab$c1)
      }
    }
  }
  
  n.contin <- ncol(tab) - 3
  if (length(latent) == 1 ){
    latent <- rep(latent, n.contin)
  } else if (length(latent) != n.contin){
    stop("Number of latent distributions does not match with the model file!")
  }
  
  B <- nrow(tab)
  
  # assign parameters to latent states:
  par.per.lat <- ifelse(latent %in% c("normal", "beta", "mises", "unif"), 2, 
                        ifelse(latent == "custom", 5, 3))
  sel.contin <- rep(1:length(latent), par.per.lat)
  X.named <- X.full <- matrix(NA, B, ncol=sum(par.per.lat))
  for(cc in 1:n.contin){
    parlist <- strsplit(tab[,3+cc], "[,()]")
    for(pp in 1:B){
      if(latent[cc] != "custom" && length(parlist[[pp]]) != par.per.lat[cc])
        stop("Number of parameters in model file (=",length(parlist[[pp]]), ")\n",
             "  does not match number of parameters of 'latent' distribution (=",par.per.lat[cc],").")
      X.named[pp,sel.contin == cc][seq_along(parlist[[pp]])] <- parlist[[pp]]
    }
  }
  
  eta.names <- sort(all.vars(parse(text=c(X.named))))
  eta.repar <- sort(na.omit(unique(c(X.named))))
  for(cc in seq_along(latent)){
    if (latent[[cc]] == "custom"){
      labels <- X.named[,which(sel.contin == cc)[1]]
      eta.names <- setdiff(eta.names, labels)
      eta.repar <- setdiff(eta.repar, labels)
    }
  }
  
  X.full <- matrix(match(X.named, eta.repar), B, sum(par.per.lat))
  X.full[is.na(X.full)] <- 0
  colnames(X.full) <- colnames(X.named) <- sel.contin
  
  if (any(X.full== 0) && all(latent != "custom"))
    stop (paste("  Check number of latent continuous parameters in model file!",
                "\n  Should be [", paste(par.per.lat,collapse=" / "),"] per line for [",
                paste(latent, collapse=" / "), "]"))
  
  # reduce duplicated states (for speed in EM algorithm)
  unique.state <- !duplicated(X.named)
  X.unique <- X.named[unique.state,,drop=FALSE] 
  X.idx <- X.full[unique.state,,drop=FALSE]   # can be recovered by get.X()
  
  # number of distinct basis distributions:
  S <- nrow(X.unique)
  map <- rep(NA, B)
  for (bb in 1:B){
    sel.state <- apply(X.unique, 1, function(xx) all(xx == X.named[bb,], na.rm = TRUE))
    map[bb] <- (1:S)[sel.state]
  }
  
  # list of S4 distributions:
  distr <- vector("list", S)
  names(distr) <- paste0("base", 1:S)
  for(s in 1:S){
    distr[[s]] <- vector("list", n.contin)
    names(distr[[s]]) <- paste0("contin", 1:n.contin)
    for(cc in 1:n.contin){
      if (latent[[cc]] == "custom"){
        idx <- as.integer(X.idx[s, sel.contin == cc])
        label <- X.unique[s,which(sel.contin == cc)[1]]
        distr[[s]][[cc]] <- make.distr(label = label, eta.idx = idx[idx != 0])
      } else {
        distr[[s]][[cc]] <- make.distr(label = latent[cc], 
                                       eta.idx = as.integer(X.idx[s, sel.contin == cc]))
      }
    }
  }
  
  list(
    # map branches to basis distributions:
    map=map, 
    # distributions:
    distr=distr, 
    # labels:
    eta.names=eta.names,
    # expression to reparameterize eta:
    eta.repar = eta.repar)
}

