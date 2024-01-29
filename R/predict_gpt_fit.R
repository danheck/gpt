#' Predictions  for GPT Models
#' 
#' Computes expected category/branch probabilities and the corresponding 
#' conditional means and quantiles/densities for the continuous latent distributions.
#' 
#' @param object GPT model fitted by \code{\link{gpt_fit}}
#' @param cat if \code{FALSE} computes expected probabilities and conditional 
#'     densities for the latent MPT branches (and not for the observable MPT categories)
#' @param dens if \code{TRUE}, returns conditional densities (instead of quantiles) for each category/branch
#' @param group select group by an index, e.g., \code{group=1} (if mutliple groups were fitted)
#' @param dim only for multivariate continuous data: dimension for prediction
#' @param quantiles which quantiles to predict
#' @param prec number of evaluations of the GPT density to compute conditional 
#'     means/quantiles for the latent distributions
#' @param ... ignored
#' 
#' @examples 
#' \dontrun{
#' # generate data
#' n <- c(targets=75, lures=75)     # number of items
#' theta <- c(do=.6,dn=.4, g=.5)          # MPT parameters
#' eta <- c(mu=400, sig=50, lambda_do=300, 
#'          lambda_go=500, lambda_gn=500, 
#'          lambda_dn=300)          # exGaussian parameters
#' file <- paste0(path.package("gpt"), "/models/2htm_exgauss.txt")
#' gen <- gpt_gen(n=n, theta=theta, eta=eta, latent="exgauss", file=file)
#' 
#' # fit GPT
#' fit <- gpt_fit(x=gen$x, y=gen$y, latent="exgauss", file=file, 
#'                restrictions=list("do=dn", "lambda_do=lambda_dn", 
#'                                  "lambda_go=lambda_gn"))
#'                                  
#' # Predictions for MPT categories:
#' predict(fit)
#' 
#' # Predictions for latent MPT branches:
#' p <- predict(fit, cat=FALSE, dens=TRUE)
#' yy <- as.numeric(colnames(p[,-(1:5)]))
#' plot(yy, p[1,-(1:5)], main="2HTM", type="l")
#' lines(yy, p[3,-(1:5)], col=2)
#' legend("topright", col=1:2, lty=1, c("Detect","Guess"))
#' }
#' 
#' @export
predict.gpt_fit <- function(object, cat = TRUE, dens = FALSE, group, dim = 1, 
                            quantiles = c(.1,.3,.5,.7,.9), prec=500, ...){
  
  object <- gpt_fit_grouped(object, group = group)
  yy <- matrix(colMeans(object$data$y), prec, ncol(object$data$y), byrow = TRUE)
  yy[,dim] <- seq(min(object$data$y[,dim]), 
                  max(object$data$y[,dim]), length.out =  prec)
  theta <- object$fit.grad$par[object$gpt@theta]
  eta <- object$fit.grad$par[object$gpt@eta]
  eta.repar <- eta.reparameterize(object$fit.grad$par[object$gpt@eta], object$gpt)
  mpt <- object$gpt@mpt
  
  ############################# OBSERVABLE CATEGORIES
  
  if (cat){
    pred <- data.frame(tree=mpt@tree.names[mpt@tree.idx],
                       cat=mpt@cat.names,
                       prob=NA,mean=NA)
    pred$prob <- mpt.cat.prob(mpt = mpt, theta = theta)
    if(!dens){
      pred <- cbind(pred, matrix(NA, nrow(pred), length(quantiles)))
      colnames(pred)[4+1:length(quantiles)] <- paste0("q",quantiles*100)
    }else{
      pred <- cbind(pred, matrix(NA, nrow(pred), prec))
      colnames(pred)[4+1:prec] <- yy[,dim]
    }
    
    for(cc in 1:nrow(pred)){
      dd <- dens(object$gpt, x = rep(cc, prec), y = yy, theta= theta, eta=eta, log = FALSE)
      pp <- sum(dd) 
      ### P(MPT cat)=pp*diff(yy)[1] # density must be scaled to y-axis!
      pred[cc,"mean"] <- sum(yy[,dim] * dd)/pp
      
      if(!dens){
        idx <- findInterval(quantiles*pp, cumsum(dd))
        qq <- (yy[idx+1,dim]+yy[idx,dim])/2
        pred[cc,4+1:length(quantiles)] <- qq
      }else{
        pred[cc,4+1:prec] <- dd
      }
    }
    
    
    ############################# OBSERVABLE CATEGORIES
  }else{
    tree.idx <- mpt@tree.idx[match(mpt@cat.names[mpt@reduce.idx], mpt@cat.names)]
    pred <- data.frame(tree=mpt@tree.names[tree.idx],
                       cat=mpt@cat.names[mpt@reduce.idx],
                       branch=1:length(mpt@reduce.idx),
                       prob=NA,mean=NA)
    pred$prob <- mpt.branch.prob(mpt, theta)
    if (!dens){
      pred <- cbind(pred, matrix(NA, nrow(pred), length(quantiles)))
      colnames(pred)[5+1:length(quantiles)] <- paste0("q",quantiles*100)
    }else{
      pred <- cbind(pred, matrix(NA, nrow(pred), prec))
      colnames(pred)[5+1:prec] <- yy[,dim]
    }
    
    for(cc in 1:nrow(pred)){
      br <- object$gpt@map[cc]
      dd <- dmultivar(y = yy, distr = object$gpt@distr[[br]],   # [dim]
                      eta = eta.repar, const = object$gpt@const, log = FALSE)
      pp <- sum(dd)
      if(!dens){
        idx <- findInterval(quantiles*pp,cumsum(dd))
        qq <- (yy[idx+1,dim]+yy[idx,dim])/2
        pred[cc,"mean"] <- sum(yy[,dim]* dd)/pp
        pred[cc,5+1:length(quantiles)] <- qq
      }else{
        pred[cc,5+1:prec] <- dd
      }
    }
  }
  
  pred
}
