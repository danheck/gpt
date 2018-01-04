

### null model:
# * categories: saturated MPT model
# * continuous: one distribution for ALL observed categories


setGeneric("model.null",
           function(model){
             standardGeneric("model.null")
           })



## gpt null for continuous:
setMethod("model.null", 
          signature(model = "gpt"), 
          function(model) {
            
            gpt.null <- model
            CC <- nrow(model@mpt@reduce)
            TT <- length(model@mpt@tree.names)
            
            # saturated MPT structure:
            gpt.null@mpt <- model.sat(model@mpt)
            
            # saturated continuous structure:
            gpt.null@map <- as.integer(rep(1, CC))
            gpt.null@theta <- names(gpt.null@mpt@theta)
            gpt.null@const <- numeric()
            
            distr <- list(list())
            eta.idx <- 0
            for(cc in seq_along(model@distr[[1]])){
              label <- model@distr[[1]][[cc]]@label
              numpar <- length(model@distr[[1]][[cc]]@lower)
              eta.idx <-  max(eta.idx) + 1:numpar
              
              # only one state/ basis distribution:
              distr[[1]][[cc]] <- make.distr(label, as.integer(eta.idx))
              
            }
            gpt.null@distr <- distr
            gpt.null@eta <- paste0("eta", 1:max(eta.idx))
            
            gpt.null
          })

