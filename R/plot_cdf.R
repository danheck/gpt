
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
cdf <- function(model, 
                ks.test = TRUE,
                dim = 1, 
                group = 1, 
                ...){
  mfrow <- par()$mfrow
  mar <- par()$mar
  
  pred <- predict(model, dens=TRUE, group=group, dim = dim) # 
  yy <- as.numeric(colnames(pred[,-(1:4)]))
  dy <- diff(yy[1:2])
  
  data <- select.data(model$data, group)
  cat.names <- model$gpt@mpt@cat.names
  J <- length(cat.names)
  
  tmp <- ceiling(sqrt(J))
  N <- length(data$x)
  miny <- min(data$y[,dim])+.1
  maxy <- max(data$y[,dim])-.1
  
  
  freq <- c(table(factor(model$data$x, labels=cat.names)))
  relfreq <- unlist(by(freq, model$gpt@mpt@tree.idx, 
                       function(x) x/sum(x)))
  names(relfreq) <- names(freq)
  N.per.tree <- c(by(freq, model$gpt@mpt@tree.idx, sum))
  names(N.per.tree) <- model$gpt@mpt@tree.names
  
  new.par <- unlist(by(rep(1,nrow(pred)), 
                       model$gpt@mpt@tree.idx, 
                       function(x) c(sum(x), rep(0, sum(x)-1))))
  
  ks <- data.frame(Category=rep(NA, J), Statistic=NA, p.value=NA)
  for(j in 1:J){
    if(new.par[j] >0) 
      par(mfrow=c(1,new.par[j]))
    y.obs <- model$data$y[model$data$x == j,dim]
    sel <- yy < max(y.obs)
    e <- ecdf(y.obs)
    pp <- cumsum(pred[j,-c(1:4), drop=TRUE])*dy
    tit <- paste0(unlist(pred[j,1:2,drop=TRUE]), collapse=": ")
    tit <- ifelse(is.null(group),tit, paste0(tit, " (group=",group,")"))
    plot(NULL, ylim=0:1, xlim = range(model$data$y[,dim]), las = 1,
         ylab="Cumulative density", xlab="Continuous Variable",
         main=tit, ...)
    ee <- e(yy)*relfreq[j]
    points(yy, ee, pch=16, col=adjustcolor("darkgray", alpha.f = .7), cex = .5)
    lines(yy[sel], pp[sel], lwd=2)
    
    # KS test: normalize to 1!
    if (ks.test){
      kstest <- ks.test(y.obs,  
                        stepfun(yy, c(0,pp/pred[j,"prob"])))
      legend("topleft", legend=paste0("D =",round(kstest$statistic,2), 
                                      "; p =",round(kstest$p.value,3)))
      ks$Category[j] <- cat.names[j]
      ks$Statistic[j] <- kstest[["statistic"]]
      ks$p.value[j] <- kstest[["p.value"]]
    }
  }
  
  par(mfrow=mfrow, mar=mar)
  if(ks.test)
    return(ks)
}
