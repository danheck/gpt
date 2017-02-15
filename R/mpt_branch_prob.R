


# standard MPT model equations: 
# theta -> branch probabilities
mpt.branch.prob <- function(mpt, theta){
  
  # constant parameters:
  tt <- theta.idx <- mpt@theta
  # free parameters:
  tt[theta.idx == -.5] <- theta
  # equality constraints:
  tt[theta.idx < -.5] <- tt[ -theta.idx[theta.idx < -.5]]

  
  B <- nrow(mpt@a)
  P1 <- length(tt)
  tmp <- matrix(tt, B, P1, byrow=T)^mpt@a * matrix(1-tt, B, P1, byrow=T)^mpt@b 
  p.branch <- apply(tmp, 1, prod)
  
  p.branch[p.branch<0] <- 0
  p.branch[p.branch>1] <- 1
  return(p.branch)
}


# theta -> category probabilities
mpt.cat.prob <- function(mpt, theta){
  
  branch.prob <- mpt.branch.prob(mpt=mpt, theta=theta)
  cat.prob <- tapply(branch.prob, mpt@reduce.idx, sum)
  
  return(cat.prob)
}