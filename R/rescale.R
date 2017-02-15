
#' Rescale bounded variable to given range
#' 
#' Uses a linear transformation to rescale variable to a specific range (e.g., [0,1] to adhere to a beta distribution). The new minimum/maximum will be equal to the lower/upper bound.
#' 
#' @param y vector of numeric values
#' @param min numeric value defining the minimum of the target range
#' @param max numeric value defining the maximum of target range
#' 
#' @export 
rescale <- function(y, min=0, max=1){
  
  if(any(is.na(y)))
    stop("Data contain missing values.")
  if(min(y) == -Inf | min==-Inf | max == Inf | min>=max | 
     length(min) !=1 | length(max) != 1)
    stop("Check values for minimum and maximum (e.g., min<max )!")
  
  # to [0,1]
  y1 <- y-min(y)
  y2 <- y1/max(y)
  
  y2*(max-min) + min
}