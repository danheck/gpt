

fit.grad <- function(gpt, x, y, 
                     starting.values, 
                     # eta.bound.list=NULL, 
                     n.fit=2, maxit=100, print = FALSE){
  
  P1 <- length(gpt@theta)
  P2 <- length(gpt@eta)
  
  eta.lower <- get.eta.lower(gpt)
  eta.upper <- get.eta.upper(gpt)

  
  # best fitting parameters:
  loglik <- rep(NA, n.fit)  
  par.mat <- matrix(NA, n.fit, P1 + P2)
  
  # res <- res.tmp <- 
  #   list(par=rep(NA, P1+P2), value = -Inf)
  
  if(print) 
    cat("\n  optim: ")
  
  for(optim.cnt in 1:n.fit){
    if(print) 
      cat( optim.cnt, "..")
    
    if(optim.cnt == 1 & !is.null(starting.values)){
      start <- starting.values
    }else{
      # test close starting values:
      start <- c(runif(P1, .2, .8),
                 random.start(eta.lower, eta.upper, starting.values[P1+1:P2]))
    }
    # scaling of eta parameters:
    eta.scale <- abs(starting.values[(P1+1:P2)[P2>0] ])
    eta.scale[eta.scale < .001] <- .001
    try({
      res <- optim(start, gpt.ll, method="L-BFGS-B", 
                       # hessian = FALSE,
                       lower=c(rep(1e-4,  P1), eta.lower), 
                       upper=c(rep(1-1e-4,P1), eta.upper),        
                       control=list(fnscale=-1, maxit=maxit, 
                                    parscale=c(rep(1,P1), eta.scale)), 
                       yy=y, xx=x, gpt=gpt)
      
      par.mat[optim.cnt,] <- res$par
      loglik[optim.cnt] <- res$val
      # if(res.tmp$value > res$value){
      #   res <- res.tmp
      # }
    }, silent=T)
  }
  
  ll.idx <- which.max(loglik)
  if(length(ll.idx) == 0){
    ll <- NA
  }else if(max(abs(outer(loglik,loglik, "-"))) > .1){
    warning("Log-likelihood differed by more than ", .1," across optim fitting runs:",
            paste(round(loglik, 2), collapse = ", "))
  }

  res <- add.vcov(gpt, par.mat[ll.idx,], x, y)
  res
}