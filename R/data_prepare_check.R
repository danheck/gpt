
#

data.check <- function(mpt, x, y, data, group=NULL){
  
  # extract from data frame:
  if(!missing(data) && !is.null(data)){
    if(!is.data.frame(data)){
      stop("data must be a data frame")
    }

    if(!is.character(x) | !is.character(y))
      stop("If a data frame 'data' is provided, x and y must define the column labels by character strings",
           "\nfor the discrete and continuous data (e.g., x='response')!")
    x <- data[,x]
    y <- data[, y, drop=FALSE]
    if(!missing(group) && !is.null(group) && length(group) == 1){
      if(!is.character(group))
        stop("If a data frame 'data' is provided, 'group' must define the column label",
             "\nfor the grouping factor (e.g., group='ID')!")
      group <- data[,group]
    }
  }
  
  x <- as.vector(x)
  y <- as.matrix(y)
  
  if(any(is.na(x) | is.na(y))){
    stop("x and y contain missing data!")
  }
  
  if(is.null(x)  || length(x) != nrow(y)){
    stop("Data x and y are not vectors of similar size!")
  }
  
  if(!is.numeric(y))
    stop("Data x and y are not numeric vectors!")
  
  if(missing(group)){
    group <- NULL
  }else{
    group <- as.vector(group)
  }
  if(!is.null(group) && length(group)!= nrow(y) )
    stop("Index variable 'group' is not a vector of the same size as x and y")
  
  
  warn <- paste("Check definition of categorical variable x:",
                "\n  either category labels as specified in the model file or",
                "\n  numbers from 1,...J for the categories:\n", paste(mpt@cat.names, collapse="\n"))
  if(!is.numeric(x)){
    x <- match(as.character(x), as.character(mpt@cat.names))
    if(sum(is.na(x) != 0)){
      stop(warn)
    }
  }else if(min(x)<1 | max(x)> length(mpt@tree.idx) | any(round(x) !=x )){
    stop(warn )
  }
  
  # file <- normalizePath(file, mustWork = TRUE)
  
  res <- list(x=x, y=y, group=group)
  return(res)
}