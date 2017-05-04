
eta.bound.par <- function(par, lower, upper, 
                          scale = .95){
  
  if(scale>=1 | scale<=.5){
    warning("scale must be in the open interval (0.5,1)")
    scale <- .95
  }
  if(any(lower >= upper)){
    whichpar <- paste0(names(lower)[lower >= upper], collapse = ", ")
    stop("lower bounds are not smaller than upper bounds: ",
         whichpar)
  }
  
  l.inf <- lower == -Inf
  u.inf <- upper == Inf
  
  # get extreme values as an anchor:
  max.par <- abs(par[ !is.na(par) & abs(par) != Inf ])
  max.bnd <- max( abs(lower[!l.inf]), 
                  abs(upper[!u.inf]))
  extreme <- min(100, max(max.bnd, max.par*10))
  
  # new, bounded ranges:
  l.bnd <- lower ; u.bnd <- upper
  # unbounded in both directions:
  l.bnd[l.inf] <-  - extreme
  u.bnd[u.inf] <-    extreme
  #l.bnd[sel.l & ! sel.u] <-  upper[sel.l &!sel.u]
  
  # 
  par <- ifelse(par <= l.bnd, 
                scale*l.bnd + (1-scale)*u.bnd,  # convex combination
                par)
  par <- ifelse(par >= u.bnd, 
                (1-scale)*l.bnd + scale*u.bnd,  # convex combination
                par)
  
  par
}