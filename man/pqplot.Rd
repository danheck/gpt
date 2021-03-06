% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_pqplot.R
\name{pqplot}
\alias{pqplot}
\title{Plot Predicted vs. Observed Quantiles}
\usage{
pqplot(
  gpt_fit,
  cumsum = FALSE,
  quantiles = seq(0.1, 0.9, 0.2),
  ncol = 4,
  dim = 1,
  group
)
}
\arguments{
\item{gpt_fit}{a fitted GPT model (see \code{\link{gpt_fit}})}

\item{cumsum}{whether to use cumulative sums of probabilities on the x axis 
(especially usefull for trees with more than 2 categories).}

\item{quantiles}{used for continuous variable}

\item{ncol}{maximum number of plots per row}

\item{dim}{only for multivariate continuous data: dimension for prediction}

\item{group}{select group by an index, e.g., \code{group=1} (if mutliple groups were fitted)}
}
\description{
Plot predicted vs. observed quantiles of continuous variable against predicted 
and observed (cumulative) category probabilities
}
