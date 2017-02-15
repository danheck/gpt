
#' Plot Predicted and Observed Densities
#' 
#' Plot predicted against observed continuous distributions
#' @param x fitted model return from \code{\link{fit.gpt}}
#' @inheritParams predict.gpt.fit
#' @param ... further arguments passt to \code{\link{hist}}
#' @importFrom graphics hist lines par plot points text
#' @export
hist.gpt.fit <- function(x, dim=1,...){
  mfrow <- par()$mfrow
  mar <- par()$mar
  
  data <- x$data
  
  cat.names <- x$gpt@mpt@cat.names
  J <- length(cat.names)
  
  tmp <- ceiling(sqrt(J))
  N <- length(data$x)
  miny <- min(data$y)-1
  maxy <- max(data$y)+1
  
  pred <- predict(x, dens=TRUE)
  yy <- as.numeric(colnames(pred[,-(1:4)]))
  scale <- diff(yy[1:2])
  
  maxdd <-max(pred[,-c(1:4)]*N*scale)
  par(mfrow=c(tmp, tmp), mar=rep(2,4))
  for(cc in 1:J){
    # sel <- pred$cat ==  x$gpt@mpt@cat.names[cc]

    dd <- pred[cc,-c(1:4)]# *scale #/pred[cc,"prob"]
    hist(x=data$y[data$x == cc], col="gray",
         breaks=seq(miny,maxy,length=min(30,floor(N/(J*2)))),# ylim=c(0, maxdd),
         main=cat.names[cc], xlab="", freq=FALSE,...)
    lines(yy, dd, col="red")
  }
  
  par(mfrow=mfrow, mar=mar)
}
