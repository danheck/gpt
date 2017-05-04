


#' Goodness of Fit of GPT Models
#' 
#' Uses a fixed numbers of histogram bins to compare GPT predictions against the empirical distribution (only for one continuous variable and group=1).
#' 
#' @param model a fitted gpt model (see \code{\link{fit.gpt}})
#' @param intervals number of intervals used for testing
#' @param breaks vector with (increasing) boundary values of intervals used for grouping
#' @seealso \code{\link{lr.test}}
#' @export
model_fit_discrete <- function (model, 
                                intervals = 10,
                                # extend = .1, 
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
                            levels = seq_along(model$gpt@mpt@cat.names),
                            labels = model$gpt@mpt@cat.names))
  # N.per.tree <- x$gpt
  miny <- min(data$y) #- extend
  maxy <- max(data$y) #+ extend
  yy <- seq(miny, maxy, intervals)
  if (!missing(breaks) && !is.null(breaks))
    yy <- breaks
  pred <- predict(yy, dens = TRUE, group=1)
  
  yy <- as.numeric(colnames(pred[,-(1:4)]))
  scale <- diff(yy[1:2])
  
  # maxdd <- max(pred[,-c(1:4)]*N*scale)
  # 
  # for(cc in 1:J){
  #   sel <- data$x == cc
  #   hh <- hist(data$y[sel], plot = FALSE,
  #              breaks=seq(miny,maxy,length=min(25,floor(N/(J*2)))))
  #   
  #   ############ Normalization
  #   dd <- unlist(pred[cc,-c(1:4)])# *scale #/pred[cc,"prob"]
  #   switch(freq,
  #          "tree" = {
  #            stop("'tree' not working")
  #            yylab <- "Density [tree-normalized]"
  #          },
  #          "cat" = {
  #            hh$counts <- hh$density
  #            dd <- dd/sum(dd)/scale #/pred[cc,"prob"]
  #            yylab <- "Density [category-normalized]"
  #          },
  #          "freq" = {
  #            hscale <- na.omit(hh$counts/hh$density)[1]
  #            dd <- hscale* dd/sum(dd)/scale #N.per.cat[cc]*scale# *scale #/pred[cc,"prob"]
  #            yylab <- "Absolute Frequency"
  #          },
  #          stop("Type of 'freq' not supported."))
  # }
  gof <- NULL
  
  return(gof)
}