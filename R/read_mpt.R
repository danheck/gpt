
###############################
# table to MPT

make.mpt <- function(tab){
  
  branch.names <- paste(tab$tree,tab$cat,sep="_" )
  
  
  # simple extractions
  cat.names <- unique(tab$cat)
  J <- length(cat.names)
  B <- nrow(tab)
  reduce.idx <- match(tab$cat,cat.names)
  names(reduce.idx) <- paste0("br", 1:B)
  tree.names <- unique(tab$tree)
  num.tree <- length(tree.names)
  tree.vec <- match(tab$tree , tree.names)
  
  
  # JxB matrix: reduce branch to category probabilities
  reduce <- matrix(0, J, B, 
                   dimnames=list( cat=cat.names, 
                                  branch=paste0("br",1:B)))
  for(rr in 1:B)
    reduce[reduce.idx[rr], rr] <- 1
  
  
  # assignment of categories to trees:
  cat.to.tree <- matrix(tree.vec, J, B, byrow=T)*reduce
  tree.idx <- apply(cat.to.tree, 1, 
                    function(xx) unique(xx[xx!=0]))
  warn <- paste0("check definition of categories and trees:",
                 "\n  some categories are assigned to more than one tree!")
  if(!is.numeric(tree.idx)){
    print(rbind(Tree=tree.names[tree.idx], Category=cat.names))
    stop(warn)
  }
  
  ####### 2.  MPT model structure: build matrices a and b
  mpt.list <- strsplit(tab$mpt, split="*", fixed=T)  
  mpt.unique <- unique(unlist(mpt.list))
  selpar <- grep("(",mpt.unique, fixed=T, invert=T)
  P1 <- length(selpar)
  theta.names <- sort(mpt.unique[selpar])
  theta.comp <- paste0("(1-",theta.names,")")
  a <- b <- matrix(0, B, P1, dimnames=list(branch.names,theta.names))
  for(bbb in 1:B){
    
    tmp <- table(mpt.list[[bbb]])
    sel_a <- theta.names %in% rownames(tmp)
    a[bbb,sel_a] <- tmp[rownames(tmp) %in% theta.names]
    sel_b <- theta.comp %in% rownames(tmp)
    b[bbb,sel_b] <- tmp[rownames(tmp) %in% theta.comp]
    
  }
  list(a=a, b=b, 
       # reduce branches to categories:
       reduce=reduce,
       reduce.idx=reduce.idx,
       # trees:
       tree.idx=tree.idx,
       # labels:
       cat.names=cat.names,
       tree.names=tree.names)
}


