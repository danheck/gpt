#' Model Selection for GPT Models
#' 
#' Currently, works only for models fitted to a single dataset
#' 
#' @param ... fitted GPT models (see \code{\link{gpt_fit}})
#' @param gpt_fits a list of fitted GPT models that are added to \code{...}
#' 
#' @seealso \code{\link{test_nested}}
#' @export
select_gpt <- function (..., gpt_fits = list()){
  
  gpt_fits <- c(list(...),  gpt_fits)
  MM <- length(gpt_fits)
  tab <- data.frame("Model" = 1:MM)
  
  ## Model Information
  tab$theta <- unlist(lapply(gpt_fits, 
                             function (x) length(x$gpt@theta)))
  tab$eta <- unlist(lapply(gpt_fits, function (x) length(x$gpt@eta)))
  tab$N <- unlist(lapply(gpt_fits, function (x) length(x$data$x)))
  
  
  ## Information Criteria
  tab$loglik <- sapply(gpt_fits, function (x) x$fit.grad$loglik)
  best <- which.max(tab$loglik)
  pars <- tab$eta + tab$theta
  tab$delta.df <- pars[best] - pars
  tab$delta.g2 <- - 2 * (tab$loglik - tab$loglik[best])
  tab$delta.pval <- pchisq(tab$delta.g2, tab$delta.df, lower.tail = FALSE)
  tab$delta.pval[best] <- NA
  tab$AIC <- sapply(gpt_fits, gpt.aic)
  tab$delta.AIC <- tab$AIC-min(tab$AIC)
  tab$wAIC <- with(tab, exp(-.5*delta.AIC)/sum(exp(-.5*delta.AIC)))
  tab$BIC <- unlist(lapply(gpt_fits, gpt.bic))
  tab$delta.BIC <- tab$BIC-min(tab$BIC)
  tab$wBIC <- with(tab, exp(-.5*delta.BIC)/sum(exp(-.5*delta.BIC)))
  
  
  if(is.null(names(gpt_fits))){
    rownames(tab) <- paste("Model", 1:MM)
  }else{
    rownames(tab) <- names(gpt_fits)
  }
  
  tab
}




########### HELPER FUNCTIONS

# @export
# gpt.G2 <- function(model1, model2){
#   warning("Ensure that models are nested!")
#   G2 <- 2*abs(model1$fit.grad$loglik - model1$fit.grad$loglik)
#   df <- abs(length(model1$gpt@eta)+length(model1$gpt@theta) - 
#               length(model2$gpt@eta)+length(model2$gpt@theta))
#   p <- pchisq(G2, df, lower.tail=F)
#   res <- c(G2=G2, df=df, prob=p)
#   return(res)
# }

gpt.aic <- function(model){
  aic <- with(model, -2*fit.grad$loglik + 
                2*(length(model$gpt@eta)+length(model$gpt@theta)))
  return(aic)
}


gpt.bic <- function(model){
  bic <- with(model, -2*fit.grad$loglik + 
                (length(model$gpt@eta)+length(model$gpt@theta))*log(length(data$x)))
  return(bic)
}
