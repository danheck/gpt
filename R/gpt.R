#' Generalized Processing Tree Models
#'
#' Fits GPT models for multivariate data (one discrete and one or more continuous responses per trial). 
#' Assumes that distribution of continuous variable(s) is a mixture distribution with the MPT core structure defining the mixture probabilities.
#'
#' The GPT structure is implemented by an S4 class \code{gpt}, which contains the MPT structure (the S4 class \code{mpt}), a vector mapping the MPT branches to the underlying continuous distributions (\code{mapvec}), a list of univariate or multivariate basis distributions (each an S4 class \code{d.uni} with information about parameter spaces etc.), the parameter labels for \code{theta} and \code{eta}, and a vector with constant values for the parameters.
#' 
#' It is advisable to first check that a GPT model file is valid using \code{\link{read_gpt}}. 
#' Next, one can either first generate some simulated data using \code{\link{gpt_gen}} 
#' or fit data using \code{\link{gpt_fit}}. The fitting algorithm first uses an EM 
#' algorithm before maximizing the full likelihood by gradient descent. 
#' Note that restrictions on the parameter space are automatically taken into account 
#' (e.g., variances must be positive). 
#' 
#' @author Daniel W. Heck, \email{heck@@uni-mannheim.de}
# @import statmod, gamlss.dist, numDeriv
#' @import stats
#' @importFrom graphics legend curve
#' @importFrom grDevices adjustcolor
#' 
#' @docType package
#' @name gpt
NULL

