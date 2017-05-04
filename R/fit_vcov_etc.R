
add.vcov <- function(gpt, par, x, y, print = FALSE){
  
  if (print) 
    cat("Estimating standard errors by Fisher information...")
  
  names(par) <- c(gpt@theta, gpt@eta)
  vcov <- hess <- SE <- NA
  ll <- gpt.ll(par=par, gpt=gpt, yy=y, xx=x)
  try ({
    hess <- hessian(gpt.ll, par, method="Richardson",
                    gpt=gpt, yy=y, xx=x)
    vcov <- -solve(hess)
    SE <- sqrt(diag(vcov))
    dimnames(vcov) <- dimnames(hess) <- 
      list(names(par),names(par))
    names(SE) <- names(par) 
  }, silent = !print)
  
  list(par=par, loglik=ll, 
       hess=hess, vcov=vcov, SE=SE)
}