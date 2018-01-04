
# evaluate expression how to reparameterize eta
eta.reparameterize <- function(eta, gpt){
  if (is.null(names(eta))){ 
    # if (length(eta) != length(gpt@eta))
    #   stop("Length of 'eta' does not match the number of parameters in GPT model.")
    names(eta) <- gpt@eta
  }
  eta.repar <- sapply(gpt@eta.repar, 
                      function(e) eval(parse(text = e), as.list(eta)))
  eta.repar
}

