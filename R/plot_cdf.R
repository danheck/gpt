
#' Plot Predicted and Empirical Cumulative Densities
#' 
#' Plot predicted against observed cumulative densities (separately per observed category)
#' @inheritParams predict.gpt.fit
#' @param model fitted GPT model (\link{fit.gpt})
#' @param ks.test whether to perform Kolmogorov-Smirnov Tests (\link{ks.test})
#' @inheritParams qplot
#' @param ... further arguments passt to \code{\link{plot}}
#' @importFrom stats ecdf ks.test
#' @export
plot_cdf <- function(model, 
                     ks.test = TRUE,
                     dim = 1, 
                     group = 1, 
                     ...){
  mfrow <- par()$mfrow
  mar <- par()$mar

  data <- select.data(model$data, group)
  cat.names <- model$gpt@mpt@cat.names
  J <- length(cat.names)
  S1 <- length(model$gpt@theta)
  cat.prob <- mpt.cat.prob(model$gpt@mpt, model$fit.grad$par[1:S1])
  
  miny <- min(data$y[,dim])+.1
  maxy <- max(data$y[,dim])-.1
  yy <- seq(miny, maxy, length.out = 101)
  
  freq <- c(table(factor(model$data$x, levels = seq_along(cat.names), labels=cat.names)))
  relfreq <- unlist(by(freq, model$gpt@mpt@tree.idx, function(x) x/sum(x)))
  names(relfreq) <- names(freq)
  N.per.tree <- c(by(freq, model$gpt@mpt@tree.idx, sum))
  names(N.per.tree) <- model$gpt@mpt@tree.names

  ks <- data.frame(Category=rep(NA, J), Statistic=NA, p.value=NA)
  new.par <- unlist(by(rep(1,J), model$gpt@mpt@tree.idx, 
                       function(x) c(sum(x), rep(0, sum(x)-1))))
  for(j in 1:J){
    if(new.par[j] >0) par(mfrow=c(1,new.par[j]))
    y.obs <- data$y[data$x == j,dim]
    e <- ecdf(y.obs)
    pp <- cdf(model$gpt, rep(j, length(yy)), matrix(yy, ncol = 1),
              model$fit.grad$par[1:S1], model$fit.grad$par[-(1:S1)], log.p = FALSE)
    tit <- paste0("category: ", cat.names[j])
    tit <- ifelse(is.null(group),tit, paste0(tit, " (group=",group,")"))
    curve(function(y) e(y)*relfreq[j], ylim=0:1, xlim = range(data$y[,dim]), las = 1,
         ylab="Cumulative density", xlab="Continuous Variable",
         col = adjustcolor("darkgray", alpha.f = .7), lwd = 1,
         main=tit, ...)
    points(yy, e(yy)*relfreq[j], pch=16, col=adjustcolor("darkgray", alpha.f = .7), cex = .5)
    lines(yy, pp, lwd=2)
    
    # KS test: normalize to 1!
    if (ks.test){
      kstest <- ks.test(y.obs, stepfun(yy, c(0,pp/cat.prob[j])))
      legend("topleft", legend=paste0("D =",round(kstest$statistic,2), 
                                      "; p =",round(kstest$p.value,3)))
      ks$Category[j] <- cat.names[j]
      ks$Statistic[j] <- kstest[["statistic"]]
      ks$p.value[j] <- kstest[["p.value"]]
    }
  }
  
  par(mfrow=mfrow, mar=mar)
  if(ks.test) return(ks)
}
