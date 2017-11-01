
add.vcov <- function(gpt, par, x, y, print = FALSE){

  names(par) <- c(gpt@theta, gpt@eta)
  hess <- NA
  try({
    hess <- hessian(gpt.ll, par, method="Richardson", gpt=gpt, yy=y, xx=x)
    dimnames(hess) <- list(names(par), names(par))
  }, silent = !print)
  
  vcov <- vcov(hess)
  for(i in 1:length(par))
    if (all(is.na(vcov)))
      vcov <- vcov(hess, omit = i)
  
  try({
    SE <- sqrt(diag(vcov))
    names(SE) <- names(par)
  }, silent = !print)
  
  ll <- gpt.ll(par=par, gpt=gpt, yy=y, xx=x)
  list("par" = par, "loglik" = ll, "hess" = hess, "vcov" = vcov, "SE" = SE)
}

vcov <- function(hess, omit = c()){
  vcov <- hess
  vcov[,] <- NA
  try ({
    if(!is.null(omit))
      vcov[- omit, - omit] <- -solve(hess[- omit, - omit])
    else
      vcov <- -solve(hess)
  }, silent = TRUE)
  vcov
}