######################## moment estimators as starting values:

guess.moments <- function(y, label, p = .9){
  start <- switch(label,
                  "normal" = c(mean = mean(y), 
                               sd = sd(y)),
                  
                  "lognormal" = {
                    shift <- min(y)*p
                    ys <- y - shift
                    c(mean = mean(log(ys)), 
                      sd = sd(log(ys)), 
                      shift = shift)}, 
                  
                  "beta" = {
                    m <- mean(y) 
                    v <- var(y)
                    a <- ifelse(v<m*(1-m), m*(m*(1-m)/v -1), 1)
                    b <- ifelse(v<m*(1-m), (1-m)*(m*(1-m)/v -1), 1)
                    c(shape1 = a, 
                      shape2 = b)},
                  
                  "gamma" = {
                    # http://research.microsoft.com/en-us/um/people/minka/papers/minka-gamma.pdf
                    # y <- rgamma(1000, 3,  15) + 10
                    shift <- max(min(y)-(1-p),0)
                    a <- .5/(log(mean(y-shift)) - mean(log(y-shift)))
                    b <- a/mean(y-shift)
                    c(shape = max(.001, a), 
                      rate = max(.001, b),
                      shift=shift)},
                  
                  "exgauss" = {
                    # heathcote (2004) BRM
                    # y <- rnorm(1000, 100,5) + rexp(100, 1/20)
                    m1 <- mean(y)
                    m2 <- var(y)
                    m3 <- sum((y-m1)^3)/length(y-1)
                    tau <- m3^(1/3)/2
                    sig <- (m2 - tau^2)^.5
                    mu <- m1 - tau
                    if (any(is.na(c(mu,sig,tau))) || any(c(mu,sig,tau) <= 0)) {
                      tau <- p*(m2^.5)
                      sig <- tau*(1-p^2)^.5
                      mu <- m1 - tau
                    }
                    # avoid too small starting values:
                    # (exGauss usually used for RTs in sec. or msec.)
                    c(mu = max(mu, .1), 
                      sigma = max(sig, .1), 
                      tau = max(tau, .1))},
                  
                  "wald" = {
                    # heathcote (2004) BRM
                    # y <- statmod:::rinvgauss(10000, 600, 400) + 320
                    # curve(statmod:::dinvgauss(x - 320, 600, 400, 200), 0, 2000)
                    shift <- p*min(y)
                    ys <- y - shift
                    m <- mean(ys)
                    c(mean = mean(ys), shape = 1 / mean(1/ys - 1/m), shift = shift)
                    },  
                  # shift range!
                  "exwald" = {
                    # heathcote (2004) BRM
                    # y <- rexwald(2000, m=10, a=5, t=10) 
                    t <- .5*sd(y)
                    m <- sqrt((mean(y)-t)/(var(y)-t^2))
                    a <- m*(mean(y)-t)
                    c(m = m, a = a, t = t)},
                  "weibull" = {
                    # y <- rweibull(1000, 2, 5) + 100
                    c(shape = .6, 
                      scale = .5, 
                      shift = min(y)*p)},
                  NA
  )
  start
}