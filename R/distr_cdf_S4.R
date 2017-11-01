

######################## DENSITY

setGeneric("cdf", 
           function(distr, x, y, theta, eta, const, log.p = FALSE) {
             standardGeneric("cdf")
           })

# univariate latent CDF:
setMethod("cdf", signature(distr = "d.uni", 
                           y="numeric", 
                           eta="numeric",
                           const = "numeric",
                           log.p = "logical"), 
          function(distr, y, eta, const, log.p = FALSE) {
            names(eta) <- names(const) <- NULL
            ee <- distr@eta.idx
            cc <- distr@eta.idx < 0
            # free parameters:
            ee[!cc] <- eta[distr@eta.idx[!cc]]
            # constants:
            ee[cc] <- const[-distr@eta.idx[cc]]
            aa <- c(list(q=y), ee, log.p=log.p)
            do.call(paste0("p", distr@dens), args = aa)
          })

#' @importFrom stats pnorm pexp pgamma

# gpt CDF
#' @importFrom matrixStats rowLogSumExps
setMethod("cdf", signature(distr = "gpt",
                           x = "numeric",
                           y = "matrix",
                           theta = "numeric",
                           eta = "numeric",
                           log.p = "logical"),
          function(distr, x, y, theta, eta, log.p=FALSE){
            N <- length(x)
            
            # MPT structure:
            branch.prob <- mpt.branch.prob(mpt=distr@mpt, theta=theta)
            lik.branch <- t(branch.prob *  t(distr@mpt@reduce))[x,,drop=FALSE]
            # browser()
            # states:
            # only for relevant rows (improved)
            if (ncol(y) > 1){
              lik.base <-  matrix(sapply(seq_along(distr@distr),
                                         # basis distribution densities:
                                         function(s) {
                                           lik.base <- rep(-Inf, N)
                                           sel.rows <- rowSums(lik.branch[,distr@map.vec == s,drop=FALSE]) != 0
                                           sel.rows[is.na(sel.rows)] <- TRUE
                                           if(any(sel.rows))
                                             lik.base[sel.rows]  <-  p.multi(y = y[sel.rows,,drop=FALSE],
                                                                             distr=distr@distr[[s]],
                                                                             eta = eta,
                                                                             const = distr@const,
                                                                             log.p = TRUE)
                                           lik.base
                                         }), N)
            } else {
              lik.base <- matrix(sapply(sapply(distr@distr, "[[", "cont1"), 
                                        cdf, y = c(y), eta=eta, const=distr@const, log.p=TRUE), N)
            }
            lik.branch <- lik.base[,distr@map.vec, drop=FALSE] + log(lik.branch)
            ll <- rowLogSumExps(lik.branch)  # = log( exp(branch 1)+...+exp(branch B) )
            ll[ll == -Inf] <- -1e100
            if (log.p){
              ll
            }else{
              exp(ll)
            }
          })
