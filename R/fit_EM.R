
#' @importFrom numDeriv hessian
#' @importFrom stats optim
fit.EM <- function(gpt, x, y, starting.values=NULL, 
                   # eta.bound.list=NULL, 
                   n.fit=5, maxit=100, tol=.001, print = FALSE){
  
  mpt <- gpt@mpt
  P1 <- length(gpt@theta)
  P2 <- length(gpt@eta)
  eta.lower <- get.eta.lower(gpt)
  eta.upper <- get.eta.upper(gpt)
  
  starting.values <- check.input.par(par = starting.values, 
                                     names = c(gpt@theta, gpt@eta))
  

  loglik <- iters <- rep(NA, n.fit)  
  par.mat <- matrix(NA, n.fit, P1 + P2)
  
  if(print) 
    cat("  EM: ")
  
  for(EM.cnt in 1:n.fit){
    if(print) 
      cat( EM.cnt, "..")
    
    if(EM.cnt == 1 && !is.null(starting.values) && length(starting.values) > 0){
      if(P1>0 && any(starting.values[1:P1]<0 | starting.values[1:P1]>1))
        stop("Check starting values: First", P1, " values must be MPT probabilites in [0,1]!")
      
      # use user-provided starting values:
      par2 <- starting.values
    }else{
      # guess starting values
      par2 <- starting.values <- guess.start(distr = gpt, y = y)
      # c(runif(P1, .2, .8),
      #           random.start(eta.lower, eta.upper, starting.values[P1+1:P2]))
    }
    par <- rep(.5,P1+P2)
    
    
    ###################################################### 
    ## EM algorithm (Dempster, Laird, & Rubin, 1977; Hu & Batchelder, 1994)
    ######################################################
    
    cnt.iter <- 0
    theta <- eta <- numeric()
    while( (sum((par-par2)^2) > tol  | cnt.iter<13) &&  # necessary due to 11 MPT-only runs!
           cnt.iter<max(15, maxit)){
      cnt.iter <- cnt.iter+1
      # old parameters:
      par <- par2 
      if(P1>0){
        theta <- par[1:P1]
      }
      if(P2 > 0){
        eta <- par[P1 + 1: P2]
      }
      
      ##########  E- step : Expectation ###############
      
      # expectation for unobserved states Z given response (standard MPT)
      tmp <- mpt.branch.prob(mpt, theta) * t(mpt@reduce)
      E.prob <- t(tmp) / colSums(tmp)
      
      if(cnt.iter>11  && P2 > 0){
        
        # probability for continuous response y for all S latent states
        S <- length(gpt@distr)
        lik.base <- matrix(0, nrow(y), S)
        for(s in 1:S){
          
          # only compute density for states with nonzero probability:
          select.rows <- rowSums(E.prob[x,][,gpt@map.vec == s,drop=FALSE]) != 0
          lik.base[select.rows,s]  <-  d.multi(y = y[select.rows,,drop=FALSE], 
                                               distr=gpt@distr[[s]], 
                                               eta = eta, 
                                               const = gpt@const, log=FALSE)
        }
        lik.branch <- lik.base[,gpt@map.vec]

        # unobserved, complete data: state indicators (each row sums up to one)
        Z.tmp <- E.prob[x,]  * lik.branch
        
        # conditional probability by normalization
        Z <-  Z.tmp /(rowSums(Z.tmp)) #E.prob[dat$x,] * E.rt.cond
      }else{
        Z <- E.prob[x,]
      }
      
      ##########  M- step: maximize parameters (separately) ###############
      
      # 1. MPT probabilities theta: analytic solution (Hu & Batchelder, 1994)
      zz <- colSums(Z)
      if(P1 > 0){
        new.theta <-  ( colSums(zz* mpt@a) / 
                           colSums( (mpt@a+mpt@b) * zz))[mpt@theta == -.5]
        new.theta <- pmax(1e-4, pmin(1-1e-4, new.theta))
        new.theta[is.na(new.theta)] <- runif(sum(is.na(new.theta)), .3,.7)
        par2[1:P1] <- new.theta
      }
      
      
      # 2. parameters eta for continuous distributions
      success <- FALSE
      if(cnt.iter>10 && P2 > 0){
        try({
          res <- optim(eta, f.complete, method="L-BFGS-B", 
                       lower=eta.lower, upper=eta.upper,
                       control=list(fnscale=-1, maxit=200,  
                                    parscale = starting.values[P1+1:P2]), 
                       y=y, gpt=gpt, Z=Z)
          success <- TRUE
        }, silent = TRUE)
        
        if(success){
          par2[P1+1:P2] <- res$par
        }else{
          par2[P1+1:P2] <- par[P1+1:P2]*runif(P2, .8, 1.2)
        }
        # print(par2)
      }else if(cnt.iter==10){
        
      }
    }
    
    # store iteration
    par.mat[EM.cnt,] <- par2
    loglik[EM.cnt] <-  gpt.ll(par=par2, gpt=gpt, xx=x, yy=y)
    iters[EM.cnt] <-  cnt.iter
  }
  
  ######## select best EM iteration as result:
  
  # cat("\nEM:",round(par.mat, 3))
  # cat("\nEM:",round(loglik, 3))
  ll.idx <- which.max(loglik)
  if(length(ll.idx) == 0){
    ll <- NA
  }else if(max(abs(outer(loglik,loglik, "-"))) > tol*10){
    warning("Log-likelihood differed by more than ", tol*10," across EM fitting runs:",
            paste(round(loglik, 2), collapse = ", "))
  }
  ll <- loglik[ll.idx]
  par.best <- par.mat[ll.idx,]
  names(par.best) <-  c(gpt@theta, gpt@eta)

  res <- list("par" = par.best, 
              "loglik" = ll, 
              "iter" = iters[ll.idx])
  res
}


