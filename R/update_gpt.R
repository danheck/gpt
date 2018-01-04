#' Re-Fit GPT Model with Additional Constraints
#' 
#' Fits nested versions of GPT models.
#' 
#' @param object a fitted GPT model (see \code{\link{gpt_fit}})
#' @param restrictions a list of additional parameter constraints that are added 
#'     to the original model
#' @param ... further arguments passed to \code{\link{gpt_fit}} for fitting the nested GPT model
#' @method update gpt_fit   
#' 
#' @seealso  \code{\link{gpt_fit}}
#' @export
update.gpt_fit <- function(object, restrictions, ...){

  # nested <- restrict_gpt(object$gpt, restrictions)
  # get_eta(nested$gpt)
  
  # gpt_fit.core(x = object$data$x, y = object$data$y, gpt = nested, ...)
}