% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_gpt_fit_cdf.R
\name{plot.gpt_fit}
\alias{plot.gpt_fit}
\title{Plot GPT Model Predictions vs. Empirical Cumulative Densities}
\usage{
\method{plot}{gpt_fit}(x, y, ks.test = FALSE, dim = 1, ...)
}
\arguments{
\item{x}{fitted model of \code{\link{gpt_fit}}}

\item{y}{ignored}

\item{ks.test}{whether to perform Kolmogorov-Smirnov Tests (\link{ks.test})}

\item{dim}{only for multivariate continuous data: dimension for prediction}

\item{...}{further arguments passt to \code{\link{plot}}}
}
\description{
Plot predicted against empirical cumulative densities, separately for each observed category.
}
