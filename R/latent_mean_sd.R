

latent.mean.sd <- function(distr, eta){
  
  switch(distr@label)
  "exwald" = {
    par <- eta[distr@eta.idx]
    m <- par[1]
    a <- par[2]
    t <- par[3]
    mean <- t + a/m
    sd <- sqrt(t^2 + a/m^3)
  }
  
  c(mean=mean, sd=sd)
}