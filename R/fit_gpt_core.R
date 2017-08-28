
#### Fitting a single data set

# x: vector of category numbers (!)
fit.gpt.core <- function( x, y, gpt, 
                             baseline = FALSE, starting.values=NULL, 
                             # eta.lower=0.01,   eta.upper=Inf, 
                             n.fit=c(6,2), maxit=c(500, 1000), EM.tol=.001, print = FALSE){
  
   
  fit.EM <- fit.EM(gpt, x, y, starting.values=starting.values, 
                   n.fit=n.fit[1], maxit=maxit[1], tol=EM.tol, print=print)
  
  if(n.fit[2]>0){
    fit.grad <- fit.grad(gpt, x, y, starting.values=fit.EM$par, 
                         n.fit=n.fit[2], maxit=maxit[2], print=print)
    
    if (!all.equal(fit.EM$par[gpt@theta], 
                   fit.grad$par[gpt@theta], tolerance=.1))
      warning ("EM and gradient estimates for MPT parameters theta differ considerably.")
    if (!all.equal(fit.EM$par[gpt@eta], 
                   fit.grad$par[gpt@eta], tolerance=.1))
      warning ("EM and gradient estimates for eta parameters differ considerably.")
    
  }else{
    fit.grad <- add.vcov(gpt, fit.EM$par, x, y)
  }
  
  
  
  # if(any(!is.null(baseline))){
  if(any(baseline)){
      if(print) 
      cat("\nFitting baseline models ...")
    
    mod.sat <- model.sat(gpt)
    fit.sat <- fit.gpt.core(x=x,y=y,  gpt=mod.sat, 
                               n.fit=n.fit, EM.tol=EM.tol, 
                               maxit=maxit, print = print, baseline=FALSE)
    df.sat <- length(fit.EM$par) - length(fit.sat$fit.EM$par) 
    G2.sat <- 2* (fit.grad$loglik - fit.sat$fit.grad$loglik )
    
    
    mod.null <- model.null(gpt)
    fit.null <- fit.gpt.core(x=x,y=y,  gpt=mod.null, 
                                n.fit=n.fit, EM.tol=EM.tol, 
                                maxit=maxit, print = print, baseline=FALSE)
    df.null <- length(fit.EM$par) - length(fit.null$fit.EM$par)
    G2.null <- 2* (fit.grad$loglik - fit.null$fit.grad$loglik )
    
    test <- list(y.per.cat = c(loglik=fit.sat$fit.grad$loglik, G2=G2.sat, df=df.sat),
                     y.null = c(loglik=fit.null$fit.grad$loglik, G2=G2.null, df=df.null))
  }else{
    test <- NULL
  }
 
  result <- list(fit.EM = fit.EM, 
                 fit.grad = fit.grad,
                 test = test)
  
  return(result)
}