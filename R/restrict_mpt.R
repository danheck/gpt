

# fix parameters
# @param restricitons either list of restrictions or location of restrictions file
# @param theta.names MPT parameter names
restrict.mpt <- function(restrictions, theta.names){
  P1 <- length(theta.names)
  # standard= -0.5 = free parameter
  theta.fixed <- rep(-0.5, length(theta.names))
  names(theta.fixed) <- theta.names
  
  # constant parameters:
  suppressWarnings(constants <- as.numeric(theta.names))
  theta.fixed[!is.na(constants)] <- constants[!is.na(constants)]
  
  if(is.null(restrictions))
    return(theta.fixed)
  
  if(!is.list(restrictions)){
    restrictions <- as.list(readLines(restrictions))
  }
  
  for(i in 1:length(restrictions)){
    
    restrictions[[i]] <- gsub(" ", "", restrictions[[i]], fixed = TRUE)
    vec <- unlist(strsplit(restrictions[[i]], "="))
    vec.num <- suppressWarnings(as.numeric(vec))
    # constant parameters
    if(length(vec) >1){
      if(any(!is.na(vec.num))){
        const <- vec.num[!is.na(vec.num)]
        if(length(const)>1)
          stop("Restrictions contain more than a single constant!")
        theta.fixed[theta.names %in% vec[is.na(vec.num)]] <- const
        # equal parameters
      }else{
        free.idx <- (1:P1)[theta.names == vec[1]]
        for(rr in 2:length(vec))
          theta.fixed[theta.names == vec[rr]] <- - free.idx
      } 
    }
  }
  
  #   res <- list(theta.fixed=theta.fixed)
  return(theta.fixed)
}

# 
# get.theta.full <- function(theta.free, theta.fixed){
#   theta.full <- theta.fixed
#   theta.full[theta.fixed == -0.5] <- theta.free 
#   sel <- - theta.fixed[theta.fixed < -.6]
#   theta.full[theta.fixed< -.6] <- theta.full[sel]
#   
#   return(theta.full)
# }