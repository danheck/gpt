#' Plot GPT Model Predictions vs. Empirical Cumulative Densities
#' 
#' Plot predicted against empirical cumulative densities, separately for each observed category.
#' 
#' @inheritParams predict.gpt_fit
#' @inheritParams pqplot
#' @inheritParams print.gpt_fit
#' @param y ignored
#' @param ks.test whether to perform Kolmogorov-Smirnov Tests (\link{ks.test})
#' @param ... further arguments passt to \code{\link{plot}}
#' 
#' @importFrom stats ecdf ks.test
#' @export
plot.gpt_fit <- function(x, y, ks.test = FALSE, dim = 1, #group = 1, 
                         ...){
  model <- x
  # mfrow <- par()$mfrow
  # mar <- par()$mar
  
  data <- select.data(model$data, group = 1)
  cat.names <- model$gpt@mpt@cat.names
  J <- length(cat.names)
  S1 <- length(model$gpt@theta)
  cat.prob <- mpt.cat.prob(model$gpt@mpt, model$fit.grad$par[1:S1])
  
  miny <- min(data$y[,dim])# +.1
  maxy <- max(data$y[,dim])# -.1
  y <- seq(miny, maxy, length.out = 101)
  
  freq <- c(table(factor(model$data$x, levels = seq_along(cat.names), labels=cat.names)))
  relfreq <- unlist(by(freq, model$gpt@mpt@tree.idx, function(x) x/sum(x)))
  names(relfreq) <- names(freq)
  N.per.tree <- c(by(freq, model$gpt@mpt@tree.idx, sum))
  names(N.per.tree) <- model$gpt@mpt@tree.names
  
  ks <- data.frame(Category=rep(NA, J), Statistic=NA, p.value=NA)
  # new.par <- unlist(by(rep(1,J), model$gpt@mpt@tree.idx, 
  #                      function(x) c(sum(x), rep(0, sum(x)-1))))
  for(j in 1:J){
    # if(new.par[j] >0) par(mfrow=c(1,new.par[j]))
    y.obs <- data$y[data$x == j,dim]
    e <- ecdf(y.obs)
    pp <- cdf(model$gpt, rep(j, length(y)), matrix(y, ncol = 1),
              model$fit.grad$par[1:S1], model$fit.grad$par[S1 + (1:S2)], log.p = FALSE)
    tit <- paste0("category: ", cat.names[j])
    # tit <- ifelse(is.null(group),tit, paste0(tit, " (group=",group,")"))
    plot(y, sapply(y, function(yy) e(yy)*relfreq[j]), type = "l", 
         ylim=0:1, xlim = range(data$y[,dim]), las = 1, 
         ylab="Cumulative density", 
         # xlab="Continuous Variable",
         col = adjustcolor("darkgray", alpha.f = .7), lwd = 2,
         main=tit, ...)
    points(y, e(y)*relfreq[j], pch=16, col="gray", cex = .8)
    lines(y, pp, lwd=2, lty = 1)
    
    # KS test: normalize to 1!
    if (ks.test){
      kstest <- ks.test(y.obs, stepfun(y, c(0,pp/cat.prob[j])))
      legend("topleft", legend=paste0("D =",round(kstest$statistic,2), 
                                      "; p =",round(kstest$p.value,3)))
      ks$Category[j] <- cat.names[j]
      ks$Statistic[j] <- kstest[["statistic"]]
      ks$p.value[j] <- kstest[["p.value"]]
    }
  }
  
  # par(mfrow=mfrow, mar=mar)
  if(ks.test) return(ks)
}
