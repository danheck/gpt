check.input.par <- function(par, names){
  
  if(is.null(par) || length(par) == 0){
    return(numeric())
  }else if(length(par) == 1){
    par <- rep(par, length(names))
  }
  
  
  # names:
  if(!is.null(names(par))){
    reorder <- match(names, names(par))
    if(min(reorder!=1) |  any(is.na(reorder)))
      stop("Names of ", substitute(par), " do not match parameters in model file!",
           "\n    Input names = ", paste(sort(names(par)), collapse=","), 
           "\n    Model names = ", paste(names, collapse=","))
    par <- par[ reorder]
  }else if(!is.null(names) & length(par) == length(names)){
    names(par) <- names
  }
  
  # lengths:
  if(length(par) != length(names))
    stop("Length of ", substitute(par), " (= ",length(par), 
         ") does not match the number of free parameters (= ",length(names), "):",
         "\n     Parameters: ", paste(names, collapse = ", "))
  
  
  par
}