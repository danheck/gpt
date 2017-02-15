# @title S4 model class for GPTs
# 
# @description Model structure for generalized processing tree (GPT) models
# @export
#' @import methods 
gpt <- setClass("gpt",  
                   representation(mpt = "mpt",          # MPT structure
                                  map.vec = "integer",   # link MPT branches to distributions
                                  distr = "list",        # list of univariate basis distributions
                                  theta = "character",   # parameter labels
                                  eta = "character",
                                  const = "numeric"),
                   validity = function(object){
                     
                     
                     if(length(unique(map.vec)) != S | any(sort.int(unique(map.vec)) != 1:S)   )
                       stop("map.vec should contain indices 1,2,3,....,S (S=number of latent continuous distributions")
                     
                     
                     
                     
                     
                   }
)

# 
setMethod(
  f="initialize",
  signature="gpt",
  definition = function(.Object, file, latent, restrictions=NULL){
    tab <- read.file.to.tab(file)
    
    mpt <- new("mpt", file, restrictions)
    gpt <- make.gpt(tab, latent=latent)
    
    .Object@mpt <- mpt
    .Object@map.vec <- as.integer(gpt$map.vec)
    .Object@theta <- colnames(mpt@a)[mpt@theta == -.5]
    
    gpt.res <- restrict.mix(gpt, restrictions)    
    .Object@distr <- gpt.res$distr
    .Object@eta <- gpt.res$eta.names
    .Object@const <- gpt.res$const
    .Object@map.vec <- as.integer(gpt.res$map.vec)
    
    return(.Object)
  }
)

# new("gpt", file, c("normal", "normal"), restrictions)
