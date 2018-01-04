#' Goodness of Fit for GPT Models
#' 
#' The continuous variable of a GPT model is categorized into discrete bins to compute Pearsons
#' X^2 between the predicted and observed bin frequencies.
#' 
#' @inheritParams print.gpt_fit
#' @param gpt_fit a fitted GPT model (see \code{\link{gpt_fit}})
#' @param breaks a list giving the breakpoints per category or a vector, in 
#'     which case the same bounds are used for each category. By default,
#'     model-implied quantiles are used for each category.
#' @param bins number of bins used to compute model-implied boundaries/quantiles.
#             or other arguments passed to \code{\link[graphics]{hist}}.
#             Can be a list of arguments if multivariate continuous data are modelled.
#' @param statistic a vector with labels of the statistic to be computed.
#' \itemize{
#'   \item \code{"dn"} the Dzhaparidze-Nikulin statistic
#'   \item \code{"pf"} the Pearson-Fisher test (refitting the model to the binned data)
#'   \item \code{"pd"} the power-divergence statistic with parameter \code{lambda}. 
#'       Since the asymptotic distribution is not chi^2, this statistic requires a 
#'       parametric or nonparametric bootstrap (not implemented).
#' }
#' @param lambda Only relevant for \code{statistic = "pf"} and \code{"pd"}. 
#'   Lambda is the parameter of the power-divergence statistic by Read & Cressie (1988).
#'   \code{lambda=1} is equivalent to Pearson's X^2, and \code{lambda=0} is
#'   equivalent to the likelihood-ratio statistic G^2.
# ' @details 
#  For the number of cells, D'Agostino & Stephens (1986; p.70) recommend 2*n^2/5 where
#  n is then number of observations for a univariate continuous random variable.
#' 
#' @template dzhaparidze1974
#' @template read1988
#' @template dagostino1986
#' @export
test_fit <- function (gpt_fit, breaks, bins = 6, statistic = "dn", lambda = 1){
  
  statistic <- sapply(statistic, match.arg, choices = c("dn", "pf", "pd"))
  names(statistic) <- statistic
  
  if(ncol(gpt_fit$data$y) > 1)
    stop("Goodness of fit not implemented for higher dimensions of y!")
  # hists <- mapply(hist, x = gpt_fit$data$y, breaks = breaks, MoreArgs = list(plot = FALSE))
  
  J <- length(gpt_fit$gpt@mpt@cat.names)  # number of categories
  S1 <- length(gpt_fit$gpt@theta)
  par <- gpt_fit$fit.grad$par
  y <- gpt_fit$data$y
  
  # model-implied quantiles (conditional on category)
  if (missing(breaks) || is.null(breaks)){
    pp <- seq(0, 1, length.out = bins + 1)[- c(1, bins + 1)]
    catprob <- mpt.cat.prob(gpt_fit$gpt@mpt, theta = par[1:S1])
    conditional_cdf <- function(y, j) 
      cdf(gpt_fit$gpt, x = rep(j, length(y)), y = matrix(y), theta = par[1:S1], 
                eta = par[-c(1:S1)], const = gpt_fit$gpt@const, log.p = FALSE) / catprob[j]
    
    diff_p <- function(q, j, cumprob) conditional_cdf(q, j) - cumprob
    
    b.list <- vector("list", J)
    names(b.list) <- gpt_fit$gpt@mpt@cat.names
    for (j in 1:J){
      for (b in 1:(bins - 1)){
        b.list[[j]][b] <- uniroot(diff_p, 
                                  lower = 2 * min(y) - max(y), 
                                  upper = 2 * max(y) - min(y), j = j, cumprob = pp[b])$root
      }
    }
    
  } else if (!is.list(breaks) & is.vector(breaks)){
    bins <- length(breaks) + 1
    b.list <- rep(list(breaks), J)
  } else if (is.list(breaks)) {
    if (length(breaks) != J)
      stop("Length of 'breaks' must be identical to number of categories.")
    b.list <- breaks
    bins <- length(breaks[[1]]) + 1
  } else {
    stop("'breaks' must be a list or vector.")
  }
  
  y.cat <- y
  for(j in 1:J){
    sel <- gpt_fit$data$x == j
    y.cat[sel] <- findInterval(y[sel], b.list[[j]]) + 1
  }
  
  N <- table(gpt_fit$gpt@mpt@tree.idx[gpt_fit$data$x])  # number of samples per tree
  prob <- prob.bins(gpt_fit$fit.grad$par, gpt_fit$gpt, b.list)
  e <- c(N[gpt_fit$gpt@mpt@tree.idx]) * prob
  o <- table(x = factor(gpt_fit$data$x, seq(gpt_fit$gpt@mpt@cat.names)), 
             y = factor(y.cat, seq(bins)))
  
  test <- lapply(paste0("test.", statistic), function(ss) 
    do.call(ss, args = list(o=o, e=e, gpt_fit = gpt_fit, b.list = b.list, lambda = lambda)))
  list("test" = data.frame(do.call("rbind", test)), 
       "observed" = o, "expected" = e, "breaks" = b.list)
}


prob.bins <- function(par, gpt, b.list){
  J <- length(gpt@mpt@cat.names)  # number of categories
  S1 <- length(gpt@theta)         # number of theta parameters
  B <- length(b.list[[1]]) + 1    # number of bins
  
  prob <- matrix(NA, J, B, dimnames = list(gpt@mpt@cat.names, paste0("bin", 1:B)))
  for (j in 1:J){
    prob[j,1:(B-1)] <- cdf(gpt, x = rep(j, B - 1), y = matrix(b.list[[j]]),
                           theta = par[1:S1], eta = par[-(1:S1)], log.p = FALSE)
  }
  prob[,B] <- mpt.cat.prob(gpt@mpt, par[1:S1])
  prob[,2:B] <- prob[,2:B]  - prob[,1:(B - 1)]
  prob
}

########## higher dimensions (not working)
# b.mat <- as.matrix(expand.grid(b.list))
# prob.list <- list()
# for(d in 1:length(b.list)){  # dimension of y
#   prob.list[[d]] <- matrix(NA, J, B, dimnames = list(gpt@mpt@cat.names, NULL))
#   for (cc in 1:J){
#     prob.list[[d]][cc,] <- cdf(gpt, x = rep(cc, B[d]), y = matrix(b.list.interior[[d]], ncol = 1),
#                                theta = par[1:S1], eta = par[-(1:S1)], log.p = FALSE)
#   }
# }

# frequency table for binned continuous (y) & discrete (x) data
count.bins <- function(x, y, J, b.list){
  B <- sapply(b.list, length) - 1  # number of bins per dependent variable
  y.cat <- data.frame(mapply(findInterval, x = data.frame(y), vec = b.list))
  for (d in 1:ncol(y.cat))
    y.cat[,d] <- factor(y.cat[,d], levels = 1:B[[d]])
  tab <- table(data.frame(x = factor(x, levels = 1:J), y.cat))
  matrix(tab, J)
}


# Power-divergence statistic
#
# o: observed frequency (book: q = probability)
# e: expected frequency (book: p = probability)
pd <- function(o, e, lambda = 1){
  pd <- switch(as.character(lambda),
               "1" = sum((o - e)^2 / e),  # X^2
               "0" = {
                 log.oe <- ifelse(o == 0, 0, log(o / e))
                 2 * sum(o*log.oe )},  # G^2
               2 * sum(o*((o/e)^lambda - 1)) / lambda / (lambda + 1))
  
  pd
}
# check continuiuty of X^2(lambda=1) and G^2 (lambda=0):
# p <- c(.1,.2,.3,.4)
# e <- p * 50
# o <- rmultinom(1, 50, p)
# curve(sapply(x, function(l) pd(o, e,l)), .9,1.1, n=10001)

# wrapper with p-value (wrong df!)
test.pd <- function(o, e, gpt_fit = NULL, b.list = NULL, lambda = 1){
  pd <- pd(o, e, lambda = lambda)
  # df <- length(o) - length(gpt_fit$gpt@mpt@tree.names) - length(gpt_fit$fit.grad$par)
  # pval <- pchisq(pd, df, lower.tail = FALSE)
  c("statistic" = pd, "df" = NA, "p.value" = NA)
}

# multinomial loglikelihood of binned frequencies:
ll.bins <- function(par, gpt, freq, b.list){
  pp <- prob.bins(par, gpt, b.list)
  nlogp <- freq * log(pp)
  nlogp[is.nan(nlogp) | nlogp == -Inf] <- -1e20
  sum(nlogp)
}

# matrix B for Rao-Robson statistic
get.B <- function (par, gpt, b.list){
  dp.dtheta <- jacobian(prob.bins, par, gpt = gpt, b.list = b.list)
  p <- c(prob.bins(par, gpt = gpt, b.list = b.list))
  1/sqrt(p) * dp.dtheta
}

################ Rao and Robson (1974): chi-square with df = M-1 
# requires expected Fihser information of sample size 1!
#' @importFrom numDeriv jacobian
test.rr <- function(o, e, gpt_fit, b.list, lambda = NULL){
  B <- get.B(gpt_fit$fit.grad$par, gpt_fit$gpt, b.list) 
  # information matrix of raw data (estimate of expected Fisher information):
  I <- - gpt_fit$fit.grad$hess / nrow(gpt_fit$data$y)  
  
  # checks whether precision decreases with binned data:
  # N <- table(gpt_fit$gpt@mpt@tree.idx[gpt_fit$data$x])  # number of samples per tree
  # J <- - hessian(ll.bins, gpt_fit$fit.grad$par, gpt = gpt_fit$gpt, freq=o, b.list = b.list)  # fisher info for binned data
  # rbind(raw = gpt_fit$fit.grad$SE,
  #       invI = 1/sqrt(diag(I*nrow(gpt_fit$data$y))),
  #       binned = sqrt(diag(solve(J))),
  #       invJ = 1/sqrt(diag(J)),
  #       invBtB = 1/sqrt(diag(t(B) %*% B * length(o) )))
  # J / length(o)  - t(B) %*% B
  
  V <- c((o - e) / sqrt(e))
  # D'Agostino & Stephens:
  Q <- diag(length(o)) + B %*% solve(I - t(B) %*% B) %*% t(B)
  # McCulloch 1985:
  # Q <- solve(diag(length(o)) + B %*% solve(- I) %*% t(B))
  # Balakrishnan, Voinov & Nikulin (2013) & D'Agostino & Stephens (Version 2):
  # rr <- pd(o, e, 1) - t(V) %*% B %*% solve(I - t(B) %*% B) %*% t(B) %*% V
  
  # Nikulin & Chimitova (2017) p. 8 (1.17) & p. 13 (1.24)
  # q <- sqrt(c(prob.bins(gpt_fit$fit.grad$par, gpt = gpt_fit$gpt, b.list = b.list)))
  # Q <- diag(length(o)) - q %*% t(q) - B %*% solve(I) %*% t(B)  # called "G" in paper
  
  rr <- c(t(V) %*% Q %*% V)
  df <- length(o) - length(gpt_fit$gpt@mpt@tree.names)
  c("statistic" = rr, "df" = df, "p.value" = pchisq(rr, df, lower.tail = FALSE))
}

################## Dzhaparidze and Nikulin (1974): df = M-p-1 
test.dn <- function (o, e, gpt_fit, b.list, lambda = NULL){
  dn <- NA
  try({
    B <- get.B(gpt_fit$fit.grad$par, gpt_fit$gpt, b.list) 
    V <- c((o - e) / sqrt(e))
    Q <- diag(length(o)) - B %*% solve(t(B) %*% B) %*% t(B)
    dn <- c(t(V) %*% Q %*% V)
  })
  df <- length(o) - length(gpt_fit$fit.grad$par) - length(gpt_fit$gpt@mpt@tree.names)
  c("statistic" = dn, "df" = df, "p.value" = pchisq(dn, df, lower.tail = FALSE))
}


################## Pearson-Fisher (re-fitting with discrete bins)
test.pf <- function(o, e, gpt_fit, b.list, lambda = lambda){
  N <- table(gpt_fit$gpt@mpt@tree.idx[gpt_fit$data$x])  # number of samples per tree
  # minimum X^2 estimation
  pd.tmp <- function(par){
    pp <- prob.bins(par, gpt = gpt_fit$gpt, b.list = b.list)
    ee <- c(N[gpt_fit$gpt@mpt@tree.idx]) * pp
    pd <- pd(o = o, e = ee, lambda = lambda)
    if (is.na(pd) || pd == Inf) pd <- 1e20
    pd
  }
  S1 <- length(gpt_fit$gpt@theta)
  oo <- optim(gpt_fit$fit.grad$par, pd.tmp,
              lower = c(rep(1e-8, S1), get.eta.lower(gpt_fit$gpt)),
              upper = c(rep(1-1e-8, S1), get.eta.upper(gpt_fit$gpt)),
              method = "L-BFGS-B", control=list(maxit = 5000))
  if (oo$convergence != 0) warning("optim did not convergence when refitting with binned data!")
  df <- length(o) - length(gpt_fit$fit.grad$par) - length(gpt_fit$gpt@mpt@tree.names)
  pval <- pchisq(oo$value, df, lower.tail = FALSE)
  c("statistic" = oo$value, "df" = df, "p.value" = pval)
}