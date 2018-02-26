#' Plot Predicted and Observed Densities
#' 
#' Plot predicted against observed continuous distributions (only for group = 1)
#' 
#' @param x fitted GPT model as returned by \code{\link{gpt_fit}}
#' @param freq how to normalize histogram and density. 
#'     Either per category (\code{freq = "cat"}), per tree (\code{freq = "tree"}), 
#'     or in absolute frequencies (\code{freq = "freq"}).
# ' @param group group to plot (not tested!)
#' @inheritParams predict.gpt_fit
#' @param ... further arguments passt to \code{\link{hist}}
#' @importFrom graphics hist lines par plot points text
#' @export
hist.gpt_fit <- function(x, dim=1, freq = "cat", # group = 1,
                     ...){
  mfrow <- par()$mfrow
  mar <- par()$mar
  
  group <- 1
  if (!is.null(x$group)){
    stop ("group not implemented")
  } else {
    data <- x$data
  }
  
  cat.names <- x$gpt@mpt@cat.names
  J <- length(cat.names)
  
  tmp <- ceiling(sqrt(J))
  N <- length(data$x)
  N.per.cat <- table(factor(data$x, 
                            levels = seq_along(x$gpt@mpt@cat.names),
                            labels = x$gpt@mpt@cat.names))
  # N.per.tree <- x$gpt
  miny <- min(data$y) #-1
  maxy <- max(data$y) #+1
  pred <- predict(x, dens = TRUE, group = 1)
  
  yy <- as.numeric(colnames(pred[,-(1:4)]))
  scale <- diff(yy[1:2])
  
  maxdd <- max(pred[,-c(1:4)]*N*scale)
  
  dots <- list(...)
  plot.args <- list(col="gray", border = "gray", xlab="", las = 1)
  if ("xlab" %in% names(dots)){
    plot.args$xlab <- dots$xlab
    dots$xlab <- NULL
  }
  # plot.args <- c(plot.args, dots)
  
  par(mfrow=c(tmp, tmp), mar=c(ifelse(plot.args$xlab == "", 1, 4),4,3,.5))
  for(cc in 1:J){
    sel <- data$x == cc
    # hh <- hist(data$y[sel], plot = FALSE, args = dots#,
               # breaks=seq(miny,maxy,length=min(25,floor(N/(J*2))))
               # )
    hh <- do.call("hist", args = c(list(x = data$y[sel], plot = FALSE), dots))
    
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
    plot.args$x = hh
    plot.args$ylab <- yylab
    plot.args$main <- cat.names[cc] 
    plot.args$ylim <- range(0, dd, hh$counts)
    do.call("plot", args = plot.args)
    # plot(hh,  col="gray", border = "gray", main=cat.names[cc], xlab=xlab, 
    #      ylim=range(0, dd, hh$counts),
    #      ylab = yylab, las = 1, ...)
    lines(yy, dd, col="black", lwd = 2)
  }
  
  par(mfrow=mfrow, mar=mar)
}
