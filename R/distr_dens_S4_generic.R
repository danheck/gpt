######################## S4 class for densities

setGeneric("dens", 
           function(distr, x, y, theta, eta, const, log=FALSE) {
             standardGeneric("dens")
           })

# univariate densities (for latent components):
setMethod("dens", signature(distr = "contin", 
                            x = "missing",
                            y = "numeric", 
                            theta = "missing",
                            eta = "numeric",
                            const = "numeric",
                            log = "logical"), 
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
                            y = "missing",
                            theta = "numeric",
                            eta = "missing",
                            const = "missing",
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
                            const = "missing",
                            log = "logical"),
          function(distr, x, y, theta, eta, log=FALSE){
            
            # browser()
            N <- length(x)
            
            # MPT structure:
            branch.prob <- mpt.branch.prob(mpt=distr@mpt, theta=theta)
            lik.branch <- t(branch.prob *  t(distr@mpt@reduce))[x,,drop=FALSE]
            
            eta.repar <- eta.reparameterize(eta, distr)
            
            # states:
            # for(ss in seq_along(distr@distr)){
            #   sel.branch <- distr@map == ss
            #   sel.rows <- rowSums(lik.branch[,sel.branch,drop=FALSE]) != 0
            #   d.contin <- dmultivar(y = y[sel.rows,,drop=FALSE],
            #                     distr=distr@distr[[ss]],
            #                     eta = eta, const = distr@const, log=FALSE)
            #   lik.branch[sel.rows, sel.branch] <- lik.branch[sel.rows, sel.branch] * d.contin
            # }
            # browser()
            # only for relevant rows (improved)
            if (ncol(y) > 1){
              lik.base <-  matrix(sapply(seq_along(distr@distr),
                                         # basis distribution densities:
                                         function(s) {
                                           lik.base <- rep(-Inf, N)
                                           sel.rows <- rowSums(lik.branch[,distr@map == s,drop=FALSE]) != 0
                                           sel.rows[is.na(sel.rows)] <- TRUE
                                           if(any(sel.rows))
                                             lik.base[sel.rows]  <-  dmultivar(y = y[sel.rows,,drop=FALSE],
                                                                               distr=distr@distr[[s]],
                                                                               eta = eta.repar,
                                                                               const = distr@const,
                                                                               log=TRUE)
                                           lik.base
                                         }), N)
            } else {
              lik.base <- matrix(sapply(sapply(distr@distr, "[[", "contin1"), 
                                        dens, y = c(y), eta=eta.repar, const=distr@const, log=TRUE), N)
              # for selected rows:
              # matrix(sapply(seq_along(distr@distr),
              #               # basis distribution densities:
              #               function(s) {
              #                 lik.base <- rep(-Inf, N)
              #                 sel.rows <- rowSums(lik.branch[,distr@map == s,drop=FALSE]) != 0
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
            #                    function(s)  dmultivar(y = y,
            #                                         distr=distr@distr[[s]],
            #                                         eta = eta,
            #                                         const = distr@const, 
            #                                         log=TRUE)), N)
            # microbenchmark::microbenchmark(#f(),f2(),
            #                                # mpt.branch.prob(mpt=distr@mpt, theta=theta),
            #                                # t(branch.prob *  t(distr@mpt@reduce))[x,,drop=FALSE], 
            #                                times = 1000)
            
            lik.branch <- lik.base[,distr@map,drop=FALSE] + log(lik.branch)
            ll <- rowLogSumExps(lik.branch)  # = log( exp(branch_1)+...+exp(branch_B) )
            
            ll[is.na(ll) | ll == -Inf] <- -1e100
            if (!log) ll <- exp(ll)
            ll
          })

