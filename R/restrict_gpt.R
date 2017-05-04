

restrict.mix <- function(gpt, restrictions){
  
  eta.names <- gpt$eta.names
  P2 <- length(eta.names)
  eta.idx.old <- eta.idx.new <- 1:P2
  
  # coding:
  # 1,2,...    = free parameter:        eta[eta.idx]
  # -1,-2,...  = constant parameter:    const[-eta.idx]
  
  
  # constant parameters (due to naming, e.g., "0.4"):
  suppressWarnings(eta.const <- as.numeric(eta.names))
  sel <- !is.na(eta.const)
  num.const <- sum(sel)
  const <- eta.const[sel]
  if(any(sel)){
    eta.idx.new[sel] <- - (1:num.const)
  }
  
  
  # check list with restrictions:
  if(!is.null(restrictions)){
    
    for(i in 1:length(restrictions)){
      
      restrictions[[i]] <- gsub(" ", "", restrictions[[i]], fixed = TRUE)
      vec <- unlist(strsplit(restrictions[[i]], "="))
      vec.num <- suppressWarnings(as.numeric(vec))
      
      if(length(vec) >1){
        # constant parameters
        if(any(!is.na(vec.num))){
          num.const <- num.const + 1
          new.const <- vec.num[!is.na(vec.num)]
          if(length(new.const)>1)
            stop("Restrictions contain more than a single constant!")
          eta.idx.new[eta.names %in% vec[is.na(vec.num)]] <- - num.const
          const[num.const] <- new.const
          # equal parameters
        }else{
          free.idx <- eta.idx.old[eta.names == vec[1]]
          for(rr in 2:length(vec))
            eta.idx.new[eta.names == vec[rr]] <-  free.idx
        }   
      }
      
    }
  }
  # CHECK:
  # restrictions; eta.names; const; eta.idx.old; eta.idx.new
  
  # re-ordering of the remaining free parameters:
  reduced.idx <- unique(eta.idx.new[eta.idx.new > 0])    # remaining free parameters (old indices)
  reduced.names <- eta.names[reduced.idx]                # remaining free parameters (labels)
  for(i in seq_along(reduced.idx)){
    eta.idx.new[eta.idx.new == reduced.idx[i]] <- i      # new indices complete!
  }
  eta.idx.new <- as.integer(eta.idx.new)
  
  # replace old by new indices in distributions:
  distr <- gpt$distr
  for(bb in seq_along(distr)){
    for(cc in seq_along(distr[[bb]])){
      distr[[bb]][[cc]]@eta.idx <- eta.idx.new[distr[[bb]][[cc]]@eta.idx]
    }
  }
  
  # eliminate identical base distributions (for speed)
  map.vec <- gpt$map.vec
  cnt <- 1
  while (max(map.vec) > cnt){
    for (cc in length(distr):(cnt+1)){
      ident <- mapply(identical, distr[[cnt]], distr[[cc]])
      if (all(ident)){
        distr[[cc]] <- NULL
        map.vec <- replace(map.vec, map.vec == cc, cnt)
        map.vec[map.vec > cc] <- map.vec[map.vec > cc] - 1
      }
    }
    cnt <- cnt + 1
  }
 names(distr) <- paste0("base", seq_along(distr)) 
  
  res <- list(distr=distr, 
              eta.names=reduced.names, 
              const=const,
              map.vec = map.vec)
  return(res)
}