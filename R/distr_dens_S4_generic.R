

######################## DENSITY

setGeneric("dens", 
           function(distr, x, y, theta, eta, const, log=FALSE) {
  standardGeneric("dens")
})

# univariate latent densities:
setMethod("dens", signature(distr = "d.uni", 
                            y="numeric", 
                            eta="numeric",
                            const = "numeric",
                            log="logical"), 
          function(distr, y, eta, const, log=FALSE) {
            names(eta) <- names(const) <- NULL
            ee <- distr@eta.idx
            cc <- distr@eta.idx < 0
            # free parameters:
            ee[!cc] <- eta[distr@eta.idx[!cc]]
            # constants:
            ee[cc] <- const[-distr@eta.idx[cc]]
            # ee <- ifelse(distr@eta.idx < 0, 
            #              const[-distr@eta.idxP1], 
            #              eta[-distr@eta.idx])
            aa <- c(list(x=y), ee, log=log)
            do.call(paste0("d", distr@dens), args = aa)
          })


#' @importFrom stats dmultinom qnorm

# MPT density (vectorized format)
# (returns matrix with responses x and latent states z)
setMethod("dens", signature(distr = "mpt",
                            x = "numeric",
                            theta = "numeric",
                            log="logical"),
          function(distr, x, theta, log=FALSE){
            
            # category probabilities:
            cat.prob <- mpt.cat.prob(mpt=distr, theta=theta)

            # observed frequencies and sample size per tree:
            freq <- as.vector(table(factor(x, levels=1:length(distr@tree.idx))))
            n <- tapply(freq, distr@tree.idx, sum)
            
            # density per tree:
            d.per.tree <- by(cbind(freq, n[distr@tree.idx], cat.prob), distr@tree.idx, 
                             function(xx) dmultinom(xx[,1], size = xx[1,2], prob = xx[,3], log = log))
            if(log){
              return(sum(d.per.tree))
            }else{
              return(prod(d.per.tree))
            }
          })


# gpt density
#' @importFrom matrixStats rowLogSumExps
setMethod("dens", signature(distr = "gpt",
                            x = "numeric",
                            y = "matrix",
                            theta = "numeric",
                            eta = "numeric",
                            log = "logical"),
          function(distr, x, y, theta, eta, log=FALSE){
            
            # browser()
            N <- length(x)

            # MPT structure:
            branch.prob <- mpt.branch.prob(mpt=distr@mpt, theta=theta)
            lik.branch <- t(branch.prob *  t(distr@mpt@reduce))[x,,drop=FALSE]
            
            # states:
            # for(ss in seq_along(distr@distr)){
            #   sel.branch <- distr@map.vec == ss
            #   sel.rows <- rowSums(lik.branch[,sel.branch,drop=FALSE]) != 0
            #   d.cont <- d.multi(y = y[sel.rows,,drop=FALSE],
            #                     distr=distr@distr[[ss]],
            #                     eta = eta, const = distr@const, log=FALSE)
            #   lik.branch[sel.rows, sel.branch] <- lik.branch[sel.rows, sel.branch] * d.cont
            # }
            # browser()
            # only for relevant rows (improved)
            if (ncol(y) > 1){
            lik.base <-  matrix(sapply(seq_along(distr@distr),
                                       # basis distribution densities:
                                       function(s) {
                                         lik.base <- rep(-Inf, N)
                                         sel.rows <- rowSums(lik.branch[,distr@map.vec == s,drop=FALSE]) != 0
                                         sel.rows[is.na(sel.rows)] <- TRUE
                                         if(any(sel.rows))
                                           lik.base[sel.rows]  <-  d.multi(y = y[sel.rows,,drop=FALSE],
                                                                           distr=distr@distr[[s]],
                                                                           eta = eta,
                                                                           const = distr@const,
                                                                           log=TRUE)
                                         lik.base
                                       }), N)
            } else {
              lik.base <- matrix(sapply(sapply(distr@distr, "[[", "cont1"), 
                                 dens, y = c(y), eta=eta, const=distr@const, log=TRUE), N)
              # for selected rows:
              # matrix(sapply(seq_along(distr@distr),
              #               # basis distribution densities:
              #               function(s) {
              #                 lik.base <- rep(-Inf, N)
              #                 sel.rows <- rowSums(lik.branch[,distr@map.vec == s,drop=FALSE]) != 0
              #                 sel.rows[is.na(sel.rows)] <- TRUE
              #                 if(any(sel.rows))
              #                   lik.base[sel.rows]  <-  dens(distr@distr[[s]][[1]],
              #                                                y = y[sel.rows,],
              #                                                eta = eta,
              #                                                const = distr@const,
              #                                                log=TRUE)
              #                 lik.base
              #               }), N)
            }
            # for all rows: slightly faster for few branches
            # lik.base <- matrix(sapply(seq_along(distr@distr), 
            #                    function(s)  d.multi(y = y,
            #                                         distr=distr@distr[[s]],
            #                                         eta = eta,
            #                                         const = distr@const, 
            #                                         log=TRUE)), N)
            # microbenchmark::microbenchmark(#f(),f2(),
            #                                # mpt.branch.prob(mpt=distr@mpt, theta=theta),
            #                                # t(branch.prob *  t(distr@mpt@reduce))[x,,drop=FALSE], 
            #                                times = 1000)
            
            lik.branch <- lik.base[,distr@map.vec,drop=FALSE] + log(lik.branch)
            ll <- sum(rowLogSumExps(lik.branch))  # = log( exp(branch 1)+...+exp(branch B) )
           

            if(ll == -Inf) 
              ll <- -1e100
            if(log){
              ll
            }else{
              exp(ll)
            }
          })
              