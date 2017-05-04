
#' Likelihood-Ratio Test
#' 
#' Test whether the increase of misfit of a nested model is significant.
#' 
#' @param model the more general model (see \code{\link{fit.gpt}})
#' @param nested the nested model (i.e., a model that differs only in additional order constraints). Note that the nestedness of models is not checked internally.
#' @param EM  whether to use the expectation-maximization fit or that of the gradient-based optimization
#' @seealso \code{\link{gpt.select}}
#' @export
lr.test <- function (model, 
                     nested, 
                     EM = TRUE){
  if (EM){
    G2 <- 2* (model$fit.EM$loglik - nested$fit.EM$loglik)
  } else {
    G2 <- 2* (model$fit.grad$loglik - nested$fit.grad$loglik)
  }
  
  S.model <- length(model$fit.EM$par)
  S.nested <- length(nested$fit.EM$par)
  df <- S.model - S.nested
  c("G2" = G2, 
    "df" = df,
    "prob" = pchisq(G2, df, lower.tail = FALSE))
}