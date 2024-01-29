
select.data <- function(data, group){
  if(!is.null(data$group)){
    if(group %in% data$group){
      sel <- data$group == group
    }else{
      sel <- data$group == unique(data$group)[group]
    }
    data <- list(x=data$x[sel],
                 y=data$y[sel,,drop=FALSE],
                 group=NULL)
  }
  data
}


gpt_fit_grouped <- function(x, group){
  if (missing(group) && !is.null(x$data$group)){
    stop("Group index must be specified to get predictions.")
  } else if (!is.null(x$data$group)){
    if(group %in% x$data$group){
      sel <- x$data$group == group
    } else if(is.numeric(group)){
      sel <- x$data$group == unique(x$data$group)[group]
    } else {
      stop ("group must either be numeric or match one of the group labels")
    }
    x <- c(x$individual.fits[[group]],  # fit.EM, fit.grad, test
           list(data = select.data(x$data, group),
                gpt=x$gpt, input = x$input))
    class(x) <- "gpt_fit"
  }
  x
}