#' Goodness of Fit for GPT Models
#' 
#' Continuous variables are categorized into discrete bins to compute Pearson's
#' X^2 between the predicted and observed bin frequencies.
#' 
#' @param model fitted GPT model (\link{fit.gpt}).
#' @param breaks a vector giving the breakpoints or other arguments passed to \code{\link[graphics]{hist}}.
#'     Can be a list of arguments if multivariate continuous data are modelled.
#' @param statistic a vector with labels of the statistic to be computed (using partial matching).
#'     \itemize{
# '     \item \code{"rr"}
#'     \item \code{"dn"} the Dzhaparidze-Nikulin statistic
#'     \item \code{"pf"} the Pearson-Fisher test (refitting the model to the binned data)
#'     \item \code{"pd"} the power-divergence statistic with parameter \code{lambda}. 
#'         Since the asymptotic distribution is not chi^2, this statistic requires a 
#'         parametric or nonparametric bootstrap (not implemented).
#'     }
#' @param lambda Only relevant for \code{statistic = "pf"} and \code{"pd"}. 
#'     Lambda is the parameter of the power-divergence statistic by Read & Cressie (1988).
#'     \code{lambda=1} is equivalent to Pearson's X^2, and \code{lambda=0} is
#'     equivalent to the likelihood-ratio statistic G^2.
#' @param bootstrap number of parametric bootstrap samples.
#' @details 
#' For the number of cells, D'Agostino & Stephens (1986; p.70) recommend 2*n^2/5 where
#' n is then number of observations for a univariate continuous random variable.
#' @template dagostino1986
# ' @template read1988
#' @export
test.gpt <- function (model, breaks = "Sturges", statistic = "dn", lambda = 1, 
                      bootstrap = 100){
  statistic <- sapply(statistic, match.arg, 
                      choices = c("dn", "pf", "pd"))
  names(statistic) <- statistic
  
  D <- ncol(model$data$y)
  if(D > 1) warning("Goodness of fit not implemented/tested for higher dimensions of y!")
  if (!is.list(breaks))
    hists <- apply(model$data$y, 2, hist, plot = FALSE, breaks = breaks)
  else
    hists <- mapply(hist, x = model$data$y, breaks = breaks, MoreArgs = list(plot = FALSE))
  
  b.list <- lapply(hists, "[[", "breaks")
  N <- table(model$gpt@mpt@tree.idx[model$data$x])  # number of samples per tree
  J <- length(model$gpt@mpt@cat.names)  # number of categories
  
  prob <- prob.bins(model$fit.grad$par, model$gpt, b.list)
  e <- c(N[model$gpt@mpt@tree.idx]) * prob
  o <- count.bins(model$data$x, model$data$y, J, b.list)
  # barplot(freq.e[c,],col=adjustcolor(2, alpha=.3))
  # barplot(freq.o[c,], add =TRUE, col=adjustcolor(4, alpha=.3))
  
  test <- sapply(paste0("test.", statistic), function(ss) 
    do.call(ss, args = list(o=o, e=e, model = model, b.list = b.list, lambda = lambda)))
  list("test" = test, "observed" = o, "expected" = e, "breaks" = b.list)
}

# debug(gpt:::prob.bins)
prob.bins <- function(par, gpt, b.list){
  J <- length(gpt@mpt@cat.names)  # number of categories
  S1 <- length(gpt@theta)  # number of theta parameters
  
  b.list.interior <- lapply(b.list, function(x) x[2:(length(x)-1)])
  b.mat <- as.matrix(expand.grid(b.list.interior))
  B <- sapply(b.list.interior, length)   # number of interior break points per dependent variable
  
  # prob.list <- list()
  # for(d in 1:length(b.list)){  # dimension of y
  #   prob.list[[d]] <- matrix(NA, J, B, dimnames = list(gpt@mpt@cat.names, NULL))
  #   for (cc in 1:J){
  #     prob.list[[d]][cc,] <- cdf(gpt, x = rep(cc, B[d]), y = matrix(b.list.interior[[d]], ncol = 1),
  #                                theta = par[1:S1], eta = par[-(1:S1)], log.p = FALSE)
  #   }
  # }
  
  cumprob <- matrix(NA, J, prod(B), dimnames = list(gpt@mpt@cat.names, NULL))
  for (cc in 1:J){
    cumprob[cc,] <- cdf(gpt, x = rep(cc, prod(B)), y = b.mat,
                        theta = par[1:S1], eta = par[-(1:S1)], log.p = FALSE)
  }
  prob.cat <- mpt.cat.prob(gpt@mpt, par[1:S1])
  prob <- cbind(cumprob[,1],
                cumprob[,2:ncol(cumprob)] - cumprob[,2:ncol(cumprob) -1],
                prob.cat - cumprob[,ncol(cumprob)])
  prob
}

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
test.pd <- function(o, e, model = NULL, b.list = NULL, lambda = 1){
  pd <- pd(o, e, lambda = lambda)
  # df <- length(o) - length(model$gpt@mpt@tree.names) - length(model$fit.grad$par)
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
test.rr <- function(o, e, model, b.list, lambda = NULL){
  B <- get.B(model$fit.grad$par, model$gpt, b.list) 
  # information matrix of raw data (estimate of expected Fisher information):
  I <- - model$fit.grad$hess / nrow(model$data$y)  
  
  # checks whether precision decreases with binned data:
  # N <- table(model$gpt@mpt@tree.idx[model$data$x])  # number of samples per tree
  # J <- - hessian(ll.bins, model$fit.grad$par, gpt = model$gpt, freq=o, b.list = b.list)  # fisher info for binned data
  # rbind(raw = model$fit.grad$SE,
  #       invI = 1/sqrt(diag(I*nrow(model$data$y))),
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
  # q <- sqrt(c(prob.bins(model$fit.grad$par, gpt = model$gpt, b.list = b.list)))
  # Q <- diag(length(o)) - q %*% t(q) - B %*% solve(I) %*% t(B)  # called "G" in paper
  
  rr <- c(t(V) %*% Q %*% V)
  df <- length(o) - length(model$gpt@mpt@tree.names)
  c("statistic" = rr, "df" = df, "p.value" = pchisq(rr, df, lower.tail = FALSE))
}

################## Dzhaparidze and Nikulin (1974): df = M-p-1 
test.dn <- function (o, e, model, b.list, lambda = NULL){
  dn <- NA
  try({
    B <- get.B(model$fit.grad$par, model$gpt, b.list) 
    V <- c((o - e) / sqrt(e))
    Q <- diag(length(o)) - B %*% solve(t(B) %*% B) %*% t(B)
    dn <- c(t(V) %*% Q %*% V)
  })
  df <- length(o) - length(model$fit.grad$par) - length(model$gpt@mpt@tree.names)
  c("statistic" = dn, "df" = df, "p.value" = pchisq(dn, df, lower.tail = FALSE))
}


################## Pearson-Fisher (re-fitting with discrete bins)
test.pf <- function(o, e, model, b.list, lambda = lambda){
  N <- table(model$gpt@mpt@tree.idx[model$data$x])  # number of samples per tree
  # minimum X^2 estimation
  pd.tmp <- function(par){
    pp <- prob.bins(par, gpt = model$gpt, b.list = b.list)
    ee <- c(N[model$gpt@mpt@tree.idx]) * pp
    pd <- pd(o = o, e = ee, lambda = lambda)
    if (is.na(pd) || pd == Inf) pd <- 1e20
    pd
  }
  S1 <- length(model$gpt@theta)
  oo <- optim(model$fit.grad$par, pd.tmp,
              lower = c(rep(1e-8, S1), get.eta.lower(model$gpt)),
              upper = c(rep(1-1e-8, S1), get.eta.upper(model$gpt)),
              method = "L-BFGS-B", control=list(maxit = 5000))
  if (oo$convergence != 0) warning("optim did not convergence when refitting with binned data!")
  df <- length(o) - length(model$fit.grad$par) - length(model$gpt@mpt@tree.names)
  pval <- pchisq(oo$value, df, lower.tail = FALSE)
  c("statistic" = oo$value, "df" = df, "p.value" = pval)
}