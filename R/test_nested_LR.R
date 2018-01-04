
#' Likelihood-Ratio Test
#' 
#' Test whether the increase of misfit of a nested model is significant.
#' 
#' @param gpt_fit the more general GPT model (see \code{\link{gpt_fit}})
#' @param nested a nested GPT model (i.e., a model that differs only in additional order constraints). 
#'     Note that the nestedness of models is not checked internally.
#' @inheritParams print.gpt_fit
#' 
#' @seealso \code{\link{select_gpt}}, \code{\link{test_fit}}
#' @export
test_nested <- function (gpt_fit, nested, EM = TRUE){
  if (EM){
    G2 <- 2* (gpt_fit$fit.EM$loglik - nested$fit.EM$loglik)
  } else {
    G2 <- 2* (gpt_fit$fit.grad$loglik - nested$fit.grad$loglik)
  }
  
  S.gpt_fit <- length(gpt_fit$fit.EM$par)
  S.nested <- length(nested$fit.EM$par)
  df <- S.gpt_fit - S.nested
  data.frame("G2" = G2, 
             "df" = df,
             "prob" = pchisq(G2, df, lower.tail = FALSE))
}