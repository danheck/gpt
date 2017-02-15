

setClass("mpt",  
         representation(a = "matrix",             # count  theta^a
                        b = "matrix",             # count (1-theta)^b
                        reduce = "matrix",        # rows: categories, cols: branches
                        reduce.idx = "integer",   # vector with category indices; length: #branches
                        tree.idx = "integer",     # tree assignment (length = categories)   # category labels
                        
                        cat.names = "character", 
                        tree.names = "character",
                        
                        theta = "numeric"      # parameter labels
         ),
         
         prototype(a = matrix(NA, 0,0),
                   b = matrix(NA, 0,0),
                   reduce = matrix(NA, 0,0),
                   reduce.idx = NA_integer_,
                   tree.idx = NA_integer_,
                   
                   cat.names = NA_character_,
                   tree.names = NA_character_,
                   theta = NA_real_),
         
         validity = function(object){
           if(any( ncol(object@reduce) != c(nrow(object@a), nrow(object@b)) ))
             stop("number of branches does not match across the arguments reduce.vec, a, and b")
           
           if(length(object@cat.names)-length(object@tree.names) < sum(object@theta == -.5))
             warning("Basic MPT model not identified: Number of parameters = ",
                     sum(object@theta == -.5),
                     "; number of free categories = ",
                     length(object@cat.names)-length(object@tree.names))
           
           
           # check whether probabilities sum up to one:
           tt <- runif(sum(object@theta == -.5))
           test <- t(mpt.cat.prob(object, tt)   * t(object@reduce) )
           sum.per.tree <- round(unlist(by(test,object@tree.idx, sum)),7)
           if(any(sum.per.tree !=1))
             warning("Check MPT definition. \n  Probabilities do not sum ob in trees: ",
                     tree.names[sum.per.tree!=1])
           TRUE
         }
)


setMethod(
  f="initialize",
  signature="mpt",
  
  definition = function(.Object, file, restrictions=NULL){
    tab <- read.file.to.tab(file)
    mpt <- make.mpt(tab)

    .Object@a <- mpt$a
    .Object@b <- mpt$b
    .Object@reduce <- mpt$reduce
    .Object@reduce.idx <- as.integer(mpt$reduce.idx)
    .Object@tree.idx <- as.integer(mpt$tree.idx)
    .Object@cat.names <- mpt$cat.names
    .Object@tree.names <- mpt$tree.names
    .Object@theta <- restrict.mpt(restrictions, colnames(mpt$a))

    return(.Object)
  }
)

