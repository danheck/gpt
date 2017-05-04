


#' Goodness of Fit of GPT Models
#' 
#' Uses a fixed numbers of histogram bins to compare model predictions against the empirical distribution.
#' 
#' @param model a fitted gpt model (see \code{\link{fit.gpt}})
#' @seealso \code{\link{lr.test}}
#' @export
model_fit_discrete <- function (model, 
                                bins = 10,
                                extend = .1, 
                                breaks){

  
  if (!is.null(model$group)){
    stop ("group not implemented")
  } else {
    data <- model$data
  }
  
  cat.names <- model$gpt@mpt@cat.names
  J <- length(cat.names)
  
  tmp <- ceiling(sqrt(J))
  N <- length(data$x)
  N.per.cat <- table(factor(data$x, 
                            levels = seq_along(x$gpt@mpt@cat.names),
                            labels = x$gpt@mpt@cat.names))
  # N.per.tree <- x$gpt
  miny <- min(data$y) - extend
  maxy <- max(data$y) + extend
  pred <- predict(x, dens = TRUE, group=group)
  
  yy <- as.numeric(colnames(pred[,-(1:4)]))
  scale <- diff(yy[1:2])
  
  maxdd <- max(pred[,-c(1:4)]*N*scale)
  
  dots <- list(...)
  plot.args <- list(col="gray", border = "gray", xlab="", las = 1)
  if ("xlab" %in% names(dots)){
    plot.args$xlab <- dots$xlab
    dots$xlab <- NULL
  }
  plot.args <- c(plot.args, dots)
  
  par(mfrow=c(tmp, tmp), mar=c(ifelse(plot.args$xlab == "", 1, 4),4,3,.5))
  for(cc in 1:J){
    sel <- data$x == cc
    hh <- hist(data$y[sel], plot = FALSE,
               breaks=seq(miny,maxy,length=min(25,floor(N/(J*2)))))
    
    ############ Normalization
    dd <- unlist(pred[cc,-c(1:4)])# *scale #/pred[cc,"prob"]
    switch(freq,
           "tree" = {
             stop("'tree' not working")
             yylab <- "Density [tree-normalized]"
           },
           "cat" = {
             hh$counts <- hh$density
             dd <- dd/sum(dd)/scale #/pred[cc,"prob"]
             yylab <- "Density [category-normalized]"
           },
           "freq" = {
             hscale <- na.omit(hh$counts/hh$density)[1]
             dd <- hscale* dd/sum(dd)/scale #N.per.cat[cc]*scale# *scale #/pred[cc,"prob"]
             yylab <- "Absolute Frequency"
           },
           stop("Type of 'freq' not supported."))
  }
  
  
  return(tab)
}