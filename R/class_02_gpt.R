# @title S4 model class for GPTs
# 
# @description Model structure for generalized processing tree (GPT) models
# @export
#' @import methods 
gpt <- setClass("gpt",  
                representation(mpt = "mpt",              # MPT structure
                               map = "integer",          # link MPT branches to distributions
                               distr = "list",           # list of univariate basis distributions
                               theta = "character",      # parameter labels
                               eta = "character",        # eta parameter labels, e.g.: "m1","m2"
                               eta.repar = "character",  # expressions used to reparameterize eta, e.g.:  "m1+m2"
                               const = "numeric"),
                validity = function(object){
                  
                  if(length(unique(map)) != S | any(sort.int(unique(map)) != 1:S)   )
                    stop("'map' must contain integers 1,2,3,....,S (S=number of latent continuous distributions")
                })

setMethod(
  f="initialize",
  signature="gpt",
  definition = function(.Object, file, latent, restrictions=NULL){
    
    
    mpt <- new("mpt", file, restrictions)
    .Object@mpt <- mpt
    .Object@theta <- colnames(mpt@a)[mpt@theta == -.5]
    
    tab <- read.file.to.tab(file)
    gpt <- make.gpt(tab, latent=latent, restrictions = restrictions)
    # gpt.res <- restrict.mix(gpt, restrictions)    
    # restrictions still apply to the eta.names parameters
    .Object@distr <- gpt$distr
    .Object@eta <- gpt$eta.names
    .Object@eta.repar <- gpt$eta.repar
    .Object@const <- numeric(0) #gpt$const
    .Object@map <- as.integer(gpt$map)
    
    .Object
  }
)
