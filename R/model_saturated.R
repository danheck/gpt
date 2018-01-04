

### saturated model:
# * categories: saturated MPT model
# * continuous: one distribution for each observed category


setGeneric("model.sat",
           function(model){
             standardGeneric("model.sat")
           })


# MPT structure
setMethod("model.sat", 
          signature(model = "mpt"), 
          function(model) {
            
            mpt.sat <- model
            
            # one branch per category:
            TT <- length(model@tree.names)
            CC <- nrow(model@reduce)
            mpt.sat@reduce.idx <- 1:CC
            mpt.sat@reduce <- diag(CC)
            rownames(mpt.sat@reduce) <- model@cat.names
            colnames(mpt.sat@reduce) <- paste0("br",1:CC)
            
            # #free parameters = #cats - #trees
            P <- CC-TT
            mpt.sat@theta <- rep(-.5, P)
            names(mpt.sat@theta) <- paste0("tt",1:P)
            
            # model matrices a and b:
            labels <- apply(cbind(mpt.sat@tree.names[mpt.sat@tree.idx], 
                                  mpt.sat@cat.names ), 1, paste, collapse="_")
            # start with zeros:
            mpt.sat@a <- mpt.sat@b <- 
              matrix(0, CC, P, dimnames = list(labels, names(mpt.sat@theta)))
            par.cnt <- 0
            for(tt in 1:TT){
              # select categories and parameters per tree:
              sel.cat <- grep(tt, mpt.sat@tree.idx)
              sel.par <- sel.cat[-length(sel.cat)] - par.cnt
              par.cnt <- max(sel.par)
              # theta^a   => diagonal = 1
              mpt.sat@a[sel.cat, sel.par] <- diag(1, length(sel.cat), length(sel.par))
              # (1-theta)^b => lower triagonal matrix = 1
              mpt.sat@b[sel.cat, sel.par] <- as.numeric(lower.tri(matrix(1, length(sel.cat), length(sel.par))))
            }

            mpt.sat
          })


## gpt saturated:
setMethod("model.sat", 
          signature(model = "gpt"), 
          function(model) {
            
            gpt.sat <- model
            CC <- nrow(model@mpt@reduce)
            TT <- length(model@mpt@tree.names)
            
            # saturated MPT structure:
            gpt.sat@mpt <- model.sat(model@mpt)
            
            # saturated continuous structure:
            gpt.sat@map <- 1:CC
            gpt.sat@theta <- names(gpt.sat@mpt@theta)
            gpt.sat@const <- numeric()
            
            distr <- list()
            eta.idx <- 0
            # as many states as categories:
            for(ss in 1:CC){
              distr[[ss]] <- list()
              
              # use first state distr[[1]] as reference:
              for(cc in seq_along(model@distr[[1]])){
                label <- model@distr[[1]][[cc]]@label
                numpar <- length(model@distr[[1]][[cc]]@lower)
                eta.idx <-  max(eta.idx) + 1:numpar
                distr[[ss]][[cc]] <- make.distr(label, as.integer(eta.idx))
              }
            }
            gpt.sat@distr <- distr
            gpt.sat@eta <- paste0("eta", 1:max(eta.idx))
            
            gpt.sat
          })
            
