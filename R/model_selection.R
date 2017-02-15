

#' Compare GPT Models
#' 
#' Currently, only functional for models fitted to single datasets
#' 
#' @param model.list list of fitted gpt-models (see \code{\link{fit.gpt}})
#' @export
gpt.select <- function(model.list){
  MM <- length(model.list)
  tab <- data.frame(Model.idx=1:MM)
  
  ## Model Information
  tab$theta <- unlist(lapply(model.list, 
                             function(x) length(x$gpt@theta)))
  tab$eta <- unlist(lapply(model.list, function(x) length(x$gpt@eta)))
  tab$N <- unlist(lapply(model.list, function(x) length(x$data$x)))
  
  
  ## Information Criteria
  tab$loglik <- unlist(lapply(model.list, 
                              function(x) x$fit$loglik))
  tab$AIC <- unlist(lapply(model.list, gpt.aic))
  tab$delta.AIC <- tab$AIC-min(tab$AIC)
  tab$wAIC <- with(tab, exp(-.5*delta.AIC)/sum(exp(-.5*delta.AIC)))
  tab$BIC <- unlist(lapply(model.list, gpt.bic))
  tab$delta.BIC <- tab$BIC-min(tab$BIC)
  tab$wBIC <- with(tab, exp(-.5*delta.BIC)/sum(exp(-.5*delta.BIC)))
  
  
  if(is.null(names(model.list))){
    rownames(tab) <- paste("Model", 1:MM)
  }else{
    rownames(tab) <- names(model.list)
  }
  
  return(tab)
}




########### HELPER FUNCTIONS

# @export
gpt.G2 <- function(model1, model2){
  warning("Ensure that models are nested!")
  G2 <- 2*abs(model1$fit$fit.grad$loglik - model1$fit$fit.grad$loglik)
  df <- abs(length(model1$gpt@eta)+length(model1$gpt@theta) - 
              length(model2$gpt@eta)+length(model2$gpt@theta))
  p <- pchisq(G2, df, lower.tail=F)
  res <- c(G2=G2, df=df, prob=p)
  return(res)
}

gpt.aic <- function(model){
  aic <- with(model, -2*fit$fit.grad$loglik + 
                2*(length(model$gpt@eta)+length(model$gpt@theta)))
  return(aic)
}


gpt.bic <- function(model){
  bic <- with(model, -2*fit$fit.grad$loglik + 
                (length(model$gpt@eta)+length(model$gpt@theta))*log(length(data$x)))
  return(bic)
}
