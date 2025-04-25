
#### Fitting a single data set

# x: vector of category numbers (!)
gpt_fit.core <- function(x, y, gpt, 
                         baseline = FALSE, starting.values=NULL, 
                         eta.lower=0.01, eta.upper=Inf,
                         n.fit=c(6,2), maxit=c(500, 1000), EM.tol=.001, print = FALSE){
  
  eta.lower <- get.eta.bound(gpt, lower = TRUE,  user.defined = eta.lower, warning = TRUE)
  eta.upper <- get.eta.bound(gpt, lower = FALSE, user.defined = eta.upper)
  
  fit.EM <- fit.EM(gpt, x, y, 
                   starting.values=starting.values, eta.lower=eta.lower, eta.upper=eta.upper,
                   n.fit=n.fit[1], maxit=maxit[1], tol=EM.tol, print=print)
  
  if (n.fit[2] > 0){
    fit.grad <- fit.grad(gpt, x, y, 
                         starting.values=fit.EM$par, eta.lower=eta.lower, eta.upper=eta.upper,
                         n.fit=n.fit[2], maxit=maxit[2], print=print)
    # try({
    if (is.character(all.equal(fit.EM$par[gpt@theta], fit.grad$par[gpt@theta], tolerance=.05)))
      warning ("EM and gradient estimates for MPT parameters theta differ considerably.")
    
    if (is.character(all.equal(fit.EM$par[gpt@eta], fit.grad$par[gpt@eta], tolerance=.05)))
      warning ("EM and gradient estimates for eta parameters differ considerably.")
    # })
  } else {
    fit.grad <- add.vcov(gpt, fit.EM$par, x, y)
  }
  
  
  
  # if(any(!is.null(baseline))){
  if (any(baseline)){
    if (print) cat("\nFitting baseline models ...")
    
    mod.sat <- model.sat(gpt)
    fit.sat <- gpt_fit.core(x=x,y=y,  gpt=mod.sat, eta.lower=eta.lower, eta.upper=eta.upper,
                            n.fit=n.fit, EM.tol=EM.tol, 
                            maxit=maxit, print = print, baseline=FALSE)
    df.sat <- length(fit.EM$par) - length(fit.sat$fit.EM$par) 
    G2.sat <- 2* (fit.grad$loglik - fit.sat$fit.grad$loglik )
    
    
    mod.null <- model.null(gpt)
    fit.null <- gpt_fit.core(x=x,y=y,  gpt=mod.null, eta.lower=eta.lower, eta.upper=eta.upper,
                             n.fit=n.fit, EM.tol=EM.tol, 
                             maxit=maxit, print = print, baseline=FALSE)
    df.null <- length(fit.EM$par) - length(fit.null$fit.EM$par)
    G2.null <- 2* (fit.grad$loglik - fit.null$fit.grad$loglik )
    
    test <- list(y.per.cat = c(loglik=fit.sat$fit.grad$loglik, G2=G2.sat, df=df.sat),
                 y.null = c(loglik=fit.null$fit.grad$loglik, G2=G2.null, df=df.null))
  } else {
    test <- NULL
  }
  
  list("fit.EM" = fit.EM, "fit.grad" = fit.grad, "test" = test)
}