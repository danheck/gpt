
#################### model equations, likelihoods etc.


# likelihood of RTs for complete data (latent states Z known)
f.complete <- function(eta, gpt, y, Z){ 
  
  eta.repar <- eta.reparameterize(eta, gpt)
  
  N <- nrow(y)
  lik.base <-  sapply(seq_along(gpt@distr),
                      # basis distribution densities:
                      function(s) {
                        lik.base <- rep(0, N)
                        # relevant rows:
                        sel.rows <- rowSums(Z[,gpt@map == s,drop=FALSE]) != 0
                        sel.rows[is.na(sel.rows)] <- TRUE
                        lik.base[sel.rows]  <-  dmultivar(y = y[sel.rows,,drop=FALSE],
                                                          distr = gpt@distr[[s]],
                                                          eta = eta.repar,
                                                          const = gpt@const, 
                                                          log=TRUE)
                        lik.base
                      })
  
  
  lik.base[lik.base == -Inf ] <- min(lik.base, -1e10, na.rm = TRUE)

  # from base distributions to branches:
  lik.branch <- lik.base[,gpt@map]
  ll <- sum(Z * lik.branch)
  if (is.na(ll) || ll == -Inf) ll <- -1e20
  ll
}


# log-likelihood of gpt model (wrapper for hessian and optim)
gpt.ll <- function (par, gpt, yy, xx){
  P1 <- length(gpt@theta)
  P2 <- length(gpt@eta)
  
  # reparameterization of eta: happens in "dens()"
  sum(dens(distr = gpt, x = xx, y = yy, 
           theta = par[seq.int(P1)], 
           eta = par[P1 + seq.int(P2)], log = TRUE))
}


