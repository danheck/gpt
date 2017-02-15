
#################### model equations, likelihoods etc.


# likelihood of RTs for complete data (latent states Z known)
f.complete <- function(eta, gpt, y, Z){ 
  
  N <- nrow(y)
  lik.base <-  sapply(seq_along(gpt@distr),
                      # basis distribution densities:
                      function(s) {
                        lik.base <- rep(0, N)
                        # relevant rows:
                        sel.rows <- rowSums(Z[,gpt@map.vec == s,drop=FALSE]) != 0
                        sel.rows[is.na(sel.rows)] <- TRUE
                        lik.base[sel.rows]  <-  d.multi(y = y[sel.rows,,drop=FALSE],
                                                        distr=gpt@distr[[s]],
                                                        eta = eta,
                                                        const = gpt@const, 
                                                        log=TRUE)
                        lik.base
                      })
  
  
  lik.base[lik.base == -Inf ] <- min(lik.base, -1e10, na.rm = TRUE)

  # from base distributions to branches:
  lik.branch <- lik.base[,gpt@map.vec]
  
  ll <- sum(Z * lik.branch)
  return(ll)
}


# log-likelihood of gpt model 
# (wrapper for hessian and optim)
gpt.ll <- function(par, gpt, yy, xx){
  P1 <- length(gpt@theta)
  dens(distr=gpt, x=xx, y=yy, 
       theta=par[1:P1], eta=par[(P1+1):length(par)], 
       log=TRUE)
}


