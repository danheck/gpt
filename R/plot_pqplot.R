#' Plot Predicted vs. Observed Quantiles
#' 
#' Plot predicted vs. observed quantiles of continuous variable against predicted 
#' and observed (cumulative) category probabilities
#' 
#' @param cumsum whether to use cumulative sums of probabilities on the x axis 
#'     (especially usefull for trees with more than 2 categories). 
#' @param quantiles used for continuous variable
#' @param ncol maximum number of plots per row
#' @inheritParams predict.gpt_fit
#' @inheritParams test_fit
#' 
#' @export
pqplot <- function(gpt_fit, cumsum = FALSE, quantiles = seq(.1,.9,.2),
                   ncol = 4, dim = 1, group){
  mfrow <- par()$mfrow
  mar <- par()$mar
  
  cat.names <- gpt_fit$gpt@mpt@cat.names
  gpt_fit <- gpt_fit_grouped(gpt_fit, group)
  x <- gpt_fit$data$x
  y <- gpt_fit$data$y[,dim]
  N <- length(x)
  freq <- c(table(factor(x, labels=cat.names)))
  
  # predictions
  pred <- obs <- predict(gpt_fit, quantiles=quantiles, dim=dim)
  sel.q <- 4+1:length(quantiles)
  
  N.per.tree <- c(by(freq, gpt_fit$gpt@mpt@tree.idx, sum))
  names(N.per.tree) <- gpt_fit$gpt@mpt@tree.names
  cnt <- 0
  for(tt in 1:length(N.per.tree)){
    sel <- which(gpt_fit$gpt@mpt@tree.idx == tt)
    cats <- obs$cat[sel]
    
    # observed probabilities and quantiles
    obs$prob[sel] <- freq[match(cats, names(freq))]/N.per.tree[tt]
    for(cc in 1:length(sel)){
      obs[sel[cc],sel.q] <- quantile(y[x == sel[cc]], quantiles)
    }
    
    if(cumsum){
      pred[sel,"prob"] <- cumsum(pred[sel,"prob"])
      obs[sel,"prob"] <- cumsum(obs[sel,"prob"])
    }
  }
  
  mini <- min(obs[,sel.q], pred[,sel.q], na.rm =TRUE)
  maxi <- max(obs[,sel.q], pred[,sel.q], na.rm = TRUE)
  
  bb <- min(ncol,  length(N.per.tree))
  aa <- ceiling( length(N.per.tree)/bb)
  par(mfrow=c(aa,bb), mar=c(2,2,2,.5))
  for(tt in 1:length(N.per.tree)){    
    sel <- which(gpt_fit$gpt@mpt@tree.idx == tt)
    cats <- cat.names[sel]
    
    plot(x=obs$prob[sel], obs[sel,5], type="p", xlim=0:1, col="blue",
         ylim=c(mini,maxi),
         xlab="", main=paste0(names(N.per.tree)[tt]," (N=",
                              N.per.tree[tt],")"), pch=8)
    text(x=obs$prob[sel], y= maxi-seq( 0,.1,length.out=length(cats))*(maxi-mini), cats)
    
    points(x=pred$prob[sel], pred[sel,5], pch=16)
    lines(x=pred$prob[sel], pred[sel,5])
    if(length(quantiles) > 1){
      for(i in 2:length(quantiles)){
        points(x=obs$prob[sel], obs[sel,4+i], pch=8,col="blue") 
        points(x=pred$prob[sel], pred[sel,4+i], pch=16)
        lines(x=pred$prob[sel], pred[sel,4+i])
      }}
  }
  
  par(mfrow=mfrow, mar=mar)
}