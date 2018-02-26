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
update.gpt_fit <- function(object, restrictions, #eta.lower, eta.upper, 
                           ...){
  
  nested <- restrict_gpt(object$gpt, restrictions)

  pars <- c(nested@theta, nested@eta)
  dots <- list(...)
  if (is.null(dots$eta.lower)) dots$eta.lower <- na.omit(object$input$eta.lower[nested@eta])
  if (is.null(dots$eta.upper)) dots$eta.upper <- na.omit(object$input$eta.upper[nested@eta])
  input <- list(x = object$data$x, y = object$data$y, gpt = nested, 
                starting.values = object$fit.EM$par[pars])
  fit <- do.call("gpt_fit.core",  c(input, dots))
  res <- c(fit, list(data = object$data, gpt=nested, 
                     input = list(file=object$input$file, latent=object$input$latent, 
                                  restrictions=restrictions)))
  class(res) <- "gpt_fit"
  res
}

restrict_gpt <- function(gpt, restrictions){
  nested <- gpt
  eta <- gpt@eta
  theta <- gpt@theta
  
  if (!missing(restrictions) && !is.null(restrictions)){
    for (i in 1:length(restrictions)){
      restrictions[[i]] <- gsub(" ", "", restrictions[[i]], fixed = TRUE)
      vec <- unlist(strsplit(restrictions[[i]], "="))
      K <- length(vec) - 1
      const <- suppressWarnings(as.numeric(vec[length(vec)]))
      if (!is.na(const)) vec <- vec[-length(vec)]
      replaced <- ifelse(is.na(const), vec[length(vec)], const)
      
      if (all(vec %in% theta)){
        stop("currently not working for 'theta' (MPT) parameters.")
      } else if (all(vec %in% eta)){
        for (k in 1:K){
          nested@eta.repar <- gsub(paste0("\\b", vec[k], "\\b"), # regex: boundary of words (not points!!)
                                   replaced, nested@eta.repar)
          nested@eta <- nested@eta[!vec[k] == nested@eta]
          for(b in seq_along(nested@distr)){
            for(d in seq_along((nested@distr[[b]]))){
              names(nested@distr[[b]][[d]]@lower) <- 
                gsub(vec[k], replaced, names(nested@distr[[b]][[d]]@lower), fixed = TRUE)
            }
          }
        }
      } else {
        stop("Restrictions must either contain only 'theta' parameters or 'eta' parameters.")
      }
    }
  }
  nested
}

# eta.names <- sort(all.vars(parse(text=c(X.named))))
# eta.repar <- sort(na.omit(unique(c(X.named))))
# for(cc in seq_along(latent)){
#   if (latent[[cc]] == "custom"){
#     labels <- X.named[,which(sel.contin == cc)[1]]
#     eta.names <- setdiff(eta.names, labels)
#     eta.repar <- setdiff(eta.repar, labels)
#   }
# }
