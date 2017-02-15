

########################################################
############## Adapted from Heathcote (2004)
############## 
############## Heathcote, A. (2004). Fitting Wald and ex-Wald distributions 
############## to response time data: An example using functions for the S-PLUS 
############## package. Behavior Research Methods, Instruments, & Computers, 36, 678–694.
########################################################


############## EX-GAUSS ############ 

dexgauss <- function(x,mu,sigma,tau, log = FALSE){
  if(!log){
    dd <- exp(((mu-x)/tau)+0.5*(sigma/tau)^2)*pnorm(((x-mu)/sigma)-(sigma/tau))/tau
    dd[is.na(dd) | x<0] <-  0 #1e-30
  }else{
    dd <- ((mu-x)/tau)+0.5*(sigma/tau)^2  + pnorm(((x-mu)/sigma)-(sigma/tau), log.p = TRUE) - log(tau)
    dd[is.na(dd) | x<0] <- -Inf # -1e10
  }
  dd
}
# log(dexg(100, 30, 4, 20))
# dexg(100, 30, 4, 20, log=TRUE)


# ex-gaussian cumulative density at x
pexgauss <- function(q,mu,sigma,tau){
  rtsu <- (q-mu)/sigma
  pnorm(rtsu) - exp((sigma^2/(2*tau^2))-((q-mu)/tau))*pnorm(rtsu - (sigma/tau))
}


# generates n random samples from the ex-gaussian (truncated normal)
rexgauss <- function(n,mu,sigma,tau){
  rr <- rep(-1, n)
  idx <- rep(TRUE, n)
  while(sum(idx)>0){
    idx <- rr<0
    rr[idx] <- rnorm(sum(idx),mean=mu,sd=sigma)	
  }
  rr + rexp(n,rate=1/tau) 
}

############## WALD / INVERSE GAUSS ############ 
dwald <- function(x,m,a,s=0, log = FALSE) {
  x <- x - s 
  suppressWarnings(
    if(!log){
      dd <- a*exp(-(a-m*x)^2/(2*x))/sqrt(2*pi*x^3)
    }else{
      dd <- log(a) -(a-m*x)^2/(2*x) - log(sqrt(2*pi*x^3))
    }
  )
  dd[is.na(dd)] <- ifelse(log, -Inf, 0)
  dd
}
# log(dwald(100, 3, 4, 2))
# dwald(100, 3, 4, 2, log=TRUE)

# Shifted Wald cumulative density with protection against numerical error
pwald <- function(q,m,a,s=0, log.p = FALSE) {
  q <- q - s
  sqrtw <- sqrt(q)
  k1 <- (m*q-a)/sqrtw
  k2 <- (m*q+a)/sqrtw
  p1 <- exp(2*a*m)
  p2 <- pnorm(-k2)
  bad <- (p1==Inf) | (p2==0); p <- p1*p2
  p[bad] <- (exp(-(k1[bad]^2)/2 - 0.94/(k2[bad]^2))/(k2[bad]*((2*pi)^.5)))
  if(!log.p){
    pp <- p + pnorm(k1)
  }else{
    pp <- log(p + pnorm(k1))
  }
  # pp[is.na(pp)] <- ifelse(log.p, -Inf, 0)
  pp
}

# Shifted Wald random function adapted from pp. 79-80, Dagpunar, J. (1988). 
# Principles of Random Variate Generation. Clarendon Press, Oxford.
rwald <- function(n,m,a,s=0) {
  if(length(n)>1) n <- length(n);
  if(length(m)>1 && length(m)!=n) m <- rep(m,length=n)
  if(length(a)>1 && length(a)!=n) lambda <- rep(a,length=n)
  y2 <- rchisq(n,1); y2onm <- y2/m; u <- runif(n)
  r1 <- (2*a + y2onm - sqrt(y2onm*(4*a+y2onm)))/(2*m)
  r2 <- (a/m)^2/r1
  ifelse(u < a/(a+m*r1), s+r1, s+r2)
}





############ EX-WALD ############ 

dexwald <- function(x,m,a,t, log = FALSE) {
  k <- (m^2 - (2/t))
  suppressWarnings(
    if (k < 0) {
      if(!log){
        dd <- exp(m*a - (a^2)/(2*x) - x*(m^2)/2)*
          rew(sqrt(-x*k/2),a/sqrt(2*x))/t
      }else{
        dd <- m*a - (a^2)/(2*x) - x*(m^2)/2 +
          log(rew(sqrt(-x*k/2),a/sqrt(2*x))) - log(t)
      }
    } else {
      k <- sqrt(k)
      if(!log){
        dd <- pwald(x,k,a)*exp(a*(m-k) - (x/t))/t
      }else{
        dd <- pwald(x,k,a, log.p = TRUE) +  a*(m-k) - (x/t) - log(t)
      }
    }
  )
  dd[is.na(dd)] <- ifelse(log, -Inf, 0)
  dd
}
# log(dexwald(100, 3, 4, 2))
# dexwald(100, 3, 4, 2, log=TRUE)

pexwald <- function(q,m,a,t) {
  pp <- pwald(q,m,a) - t*dexwald(q,m,a,t)
  pp[is.na(pp)] <- ifelse(log, -Inf, 0)
  pp
}
rexwald <- function(n,m,a,t) {
  rwald(n,m,a) + rexp(n,1/t)
}








############################# AUXILIARY FUNCTIONS ######################################

# Series approximation to the real (u) and imaginary (v) parts of complex 
# error function, erf(x + iy). Approximation can fail if firstblock > 20 
# but arguments allow forcing. Blocks loop will only process up to  
# maxseries terms.
uandv <- function(x,y,firstblock=20,block=0,tol=.Machine$double.eps^(2/3),maxseries=20)
{
  twoxy <- 2*x*y; xsq <- x^2; iexpxsqpi <- 1/(pi*exp(xsq))
  sin2xy <- sin(twoxy); cos2xy <- cos(twoxy)
  nmat <- matrix(rep((1:firstblock),each=length(x)),nrow=length(x))
  nsqmat <- nmat^2; ny <- nmat*y; twoxcoshny <- 2*x*cosh(ny)
  nsinhny <- nmat*sinh(ny); nsqfrac <- (exp(-nsqmat/4)/(nsqmat + 4*xsq))
  u <- (2*pnorm(x*sqrt(2))-1)+iexpxsqpi*(((1-cos2xy)/(2*x))+2*
                                           ((nsqfrac*(2*x-twoxcoshny*cos2xy+nsinhny*sin2xy))%*%rep(1,firstblock)))
  v <-iexpxsqpi*((sin2xy/(2*x))+2*((nsqfrac*(twoxcoshny*sin2xy+ 
                                               nsinhny*cos2xy))%*%rep(1,firstblock)))
  n <- firstblock;	converged <- rep(F,length(x))
  repeat {
    if ((block < 1) || (n >= maxseries)) break
    else {
      if ((n + block) > maxseries) block <- (maxseries - n)
      nmat <-matrix(rep((n+1):(n+block),each=sum(!converged)), 
                    nrow=sum(!converged))
      nsq <- nmat^2; ny <- nmat*y[!converged];
      twoxcoshny <- 2*x[!converged]*cosh(ny); nsinhny <- nmat*sinh(ny)
      nsqfrac <- (exp(-nsq/4)/(nsq + 4*xsq[!converged]))
      du <- iexpxsqpi[!converged]*((2*nsqfrac*(2*x[!converged]- 
                                                 twoxcoshny*cos2xy[!converged]+nsinhny*sin2xy[!converged]))  
                                   %*%rep(1,block))
      dv <-iexpxsqpi[!converged]*((2*nsqfrac*(twoxcoshny*sin2xy[!converged]+ 
                                                nsinhny*cos2xy[!converged]))%*%rep(1,block))
      u[!converged] <- u[!converged] + du;
      v[!converged] <- v[!converged] + dv
      converged[!converged] <- ((du < tol) & (dv < tol))
      if (all(converged)) break
    }
  }
  cbind(u,v)
}

# real part of w(z) = exp(-z^2)[1-erf(-iz)]
rew <- function(x,y,...) {
  uv <- uandv(y,x,...)
  exp(y^2-x^2)*(cos(2*x*y)*(1-uv[,1])+sin(2*x*y)*uv[,2])
}

