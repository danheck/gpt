

######################## RANDOM NUMBER GENERATION

setGeneric("rand", function(distr, n, eta=NULL, theta=NULL, const=NULL) {
  standardGeneric("rand")
})


# latent continuous basis distributions
# (returns matrix y)
setMethod("rand", signature(distr = "d.uni", 
                            n = "numeric", 
                            eta = "numeric",
                            const = "numeric"), 
          function(distr, n, eta, const){
            names(eta) <- names(const) <- NULL
            ee <- distr@eta.idx
            cc <- distr@eta.idx < 0
            # free parameters:
            ee[!cc] <- eta[distr@eta.idx[!cc]]
            # constants:
            ee[cc] <- const[-distr@eta.idx[cc]]
            if(any(ee < distr@lower, ee>distr@upper)){
              warning("Parameter values out of range!")
              rep(NaN, n)
            }else{
              aa <- c(list(n=n), ee)
              do.call(paste0("r", distr@dens), args = aa)
            }
          })

# MPT data generation 
# (returns matrix with responses x and latent states z)
setMethod("rand", signature(distr = "mpt",
                            n = "numeric",
                            theta = "numeric"),
          function(distr, n, theta){
            mpt <- distr
            
            if (any(theta<0, theta>1)){
              stop("MPT parameters theta must be in the interval [0,1].")
              rep(NaN, sum(n))
            }
            if (length(n) != length(mpt@tree.names))
              stop ("Vector 'n' with number of observations does not match number of MPT trees!")
            if (any(n<0))
              stop ("Vector 'n' with number of observations must be nonnegative!")
            
            # model details
            B <- ncol(mpt@reduce) # branches
            J <- nrow(mpt@reduce) # categories
            num.tree <- max(mpt@tree.idx)
            
            # vectors of discrete responses and latent indicators
            x <- z <- rep(NA, sum(n))
            
            # sample latent states (MPT branches)
            randnum <- runif(sum(n))
            prob.branch <- mpt.branch.prob(mpt, theta)
            for(tt in seq_along(mpt@tree.names)){
              
              if (n[tt] > 0){
                # select responses, categories, branches:
                sel.resp <- 1:n[tt] + ifelse(tt <= 1, 0, sum(n[1:(tt-1)]))
                sel.cat <- mpt@tree.idx == tt
                sel.branch <- (1:B)[mpt@reduce.idx %in% (1:J)[sel.cat]]
                
                prob.branch.tmp <- prob.branch[sel.branch]
                
                # if(round(sum(cat.prob.tmp),10) != 1)
                #   warning(paste("Check model definition! Probabilities do not add up to one within MPT tree", t))
                
                # generate latent state:
                idx <- findInterval(randnum[sel.resp], cumsum(prob.branch.tmp) ) +1
                z[sel.resp] <- sel.branch[idx]
                x[sel.resp] <- mpt@reduce.idx[z[sel.resp]]
              }
            }
            cbind(x, z)
          })


# gpt
# returns data frame with columns: 
# tree/discrete x/continuous y1...ym/state z
setMethod("rand", signature(distr = "gpt",
                            n = "numeric",
                            eta = "numeric",
                            theta = "numeric"),
          function(distr, n, eta, theta){
            gpt <- distr
            
            # MPT data generation (S4 method):
            mat <- rand(distr=gpt@mpt, n=n, theta=theta)
            x <- mat[,"x"]
            z <- mat[,"z"]
            
            # number of basis distributions:
            S <- length(gpt@distr)    
            num.tree <- max(gpt@mpt@tree.idx)
            
            y <- matrix(NA, sum(n), length(gpt@distr[[1]]))
            
            # get states:
            states <- gpt@map[z]
            
            # reparameterize eta
            eta.repar <- eta.reparameterize(eta, gpt)
            
            # sample continuous values based on latent state z
            y.list <- tapply(states, list(factor(states, levels=1:S)), 
                             function(ss){
                               if(!is.null(ss))
                                 r.multi(length(ss), gpt@distr[[ss[1]]], eta.repar, gpt@const)
                             })
            # correct assignment:
            for(ss in 1:S){
              if(sum(states == ss) > 0)
                y[states == ss,] <- y.list[[ss]]
            }
            
            data.frame(tree=gpt@mpt@tree.names[rep(1:num.tree, n)], 
                       x=gpt@mpt@cat.names[x], 
                       y=y, 
                       state=z)
          })