# guess starting values:
# fit simple distribution to pooled data
setGeneric("guess.start",
           function(distr, y){
             standardGeneric("guess.start")
           })


# guess starting values for a single univariate distribution:
setMethod("guess.start", signature(distr = "contin", 
                                   y = "numeric"), 
          function(distr, y) {
            qq <- quantile(y, c(.05, .95))
            y <- y[y>qq[1] & y<qq[2]]
            
            guess.moments <- guess.moments(y, distr@label)
              
            if(any(is.na(guess.moments))){
              guess.moments <- random.start(distr@lower, distr@upper)
            }
            guess.moments <- eta.bound.par(guess.moments, 
                                           distr@lower, distr@upper)
            
            # new univariate distribution without restrictions:
            dd <- make.distr(distr@label, 1:length(distr@eta.idx))

            # loglikelihood:
            ll <- function(par) {
              suppressWarnings(
                ll <- sum(dens(distr=dd, y=y, eta = par, const=NA_real_, log=TRUE))
              )
              if(is.na(ll) || ll == -Inf) ll <- -1e20
              ll
            }

            # fit parameters:
            oo <- optim(guess.moments, fn = ll, 
                        control = list(fnscale=-1, parscale=guess.moments), 
                        method = "L-BFGS-B",
                        lower = distr@lower, upper = distr@upper)
            
            
            # move starting points to the inner area of 
            # the allowed parameter space:
            guess <- eta.bound.par(oo$par, distr@lower, distr@upper)

            guess   
          })



# use one set of parameters for each continuous variable y:
setMethod("guess.start", signature(distr = "gpt", 
                                   y = "matrix"), 
          function(distr, y){
            
            gpt <- distr
            guess <- rep(NA, length(gpt@eta))
            
            # find starting values for each continuous variable
            for(cc in 1:ncol(y)){
              guess.cc <- guess.start(distr = gpt@distr[[1]][[cc]], 
                                      y = y[,cc,drop=TRUE])
              
              # assign fitted values to free parameters:
              for(gg in seq_along(guess.cc)){
                # basis distributions:
                for(ss in seq_along(gpt@distr)){
                  
                  # get eta index:
                  eta.idx.gg <- gpt@distr[[ss]][[cc]]@eta.idx[gg]
                  # if guess.cc[gg] is free parameter:
                  if(eta.idx.gg > 0){
                    
                    # replace free parameters by appropriate parameters (mean etc.):
                    guess[eta.idx.gg] <- guess.cc[gg]
                  }
                }
              }
            }
            
            guess <- c(runif(length(gpt@theta), .2,.8), guess)
            
            names(guess) <- c(gpt@theta, gpt@eta)
            guess
          })






######################## random, bounded starting values:

random.start <- function(eta.lower, eta.upper, starting.values){

  if(missing(starting.values) || is.null(starting.values)){
    eta.start <- runif(length(eta.lower), 
                       min = pmax(eta.lower, pmin(0.1, eta.upper-.01)), # avoid -Inf
                       max = pmin(eta.upper, pmax(100, eta.lower+.01)))  # avoid +Inf
  }else{
    # educated guess:
    guess.l <- ifelse(abs(starting.values) >1 , 
                      starting.values - 3 * abs(starting.values),
                      -10)
    guess.u <- ifelse(abs(starting.values) >1 , 
                      starting.values + 3 * abs(starting.values),
                      10)
    eta.start <- runif(length(eta.lower), 
                       min = pmax(eta.lower, pmin(guess.l, eta.upper-.01)),
                       max = pmin(eta.upper, pmax(guess.u, eta.lower+.01)))
  }
  
  eta.start
}