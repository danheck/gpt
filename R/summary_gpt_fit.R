#' Summary for Fitted GPT Models
#' 
#' Shows parameter estimates and goodness of fit.
#' 
#' @inheritParams update.gpt_fit
#' @inheritParams print.gpt_fit
#' @param ... further arguments passed to \code{\link{test_fit}}
#' 
#' @method summary gpt_fit
#' @export
summary.gpt_fit <- function(object, EM = FALSE, ...){
  cat("## Parameter estimates (", ifelse(EM, "EM", "optim"),")\n", sep = "")
  print(object, EM = EM)
  
  test <- test_fit(object, ...)
  cat("\n## Goodness of Fit (",ncol(test$expected)," bins):\n", sep = "")
  print(test$test)  
}

