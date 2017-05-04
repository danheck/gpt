
#' Check GPT Model File
#' 
#' Check model file and parameter restrictions
#' 
#' @param file path to model file
#' @param latent type of latent continuous distribution(s) (e.g., \code{"normal"}, \code{"gamma"}, \code{"exgauss"}, etc.). See \code{\link{fit.gpt}} for details 
#' @param restrictions list of restrictions for MPT and continuous parameters, e.g., \code{list("g=.5", "do=dn", "mu1=mu2=mu3")}
#' @details The model file should have one line per MPT branch, with entries separated by semicolons in the following order:
#' 
#' \code{Tree; Category; MPT equations;  Continuous Parameters 1;  Cont. Parameters 2}
#' 
#' In the model file, the latent distributions are specified by parameters corresponding to the distribution of each processing branch (separated by commas). The definition of the parameters is listed in \code{\link{fit.gpt}}. As an example for \code{latent=c("normal","normal")}:
#' 
#' \tabular{lllll}{
#' \code{#tree} \tab \code{; cat} \tab \code{; mpt} \tab \code{; normal 1}\cr\tab \code{; normal 2}\cr
#' \code{tree1} \tab \code{; hit} \tab \code{; do}           \tab \code{; m1_d,s1} \tab \code{; m2_d,s2}\cr
#' \code{tree1} \tab \code{; hit} \tab \code{; (1-do)*g}     \tab \code{; m1_g,s1} \tab \code{; m2_g,s2}\cr
#' \code{tree1} \tab \code{; miss}\tab \code{; (1-do)*(1-g)} \tab \code{; m1_g,s1} \tab \code{; m2_g,s2}\cr
#' \code{tree2} \tab \code{; fa}  \tab \code{; (1-dn)*g}     \tab \code{; m1_g,s1} \tab \code{; m2_g,s2}\cr
#' \code{tree2} \tab \code{; cr}  \tab \code{; (1-dn)*(1-g)} \tab \code{; m1_g,s1} \tab \code{; m2_g,s2}\cr
#' \code{tree2} \tab \code{; cr}  \tab \code{; dn}           \tab \code{; m1_d,s1} \tab \code{; m2_d,s2} 
#' }
#' @examples
#' file <- paste0(path.package("gpt"), "/models/2htm_2normal.txt")
#' check.gpt(file=file, latent=c("normal", "normal"))
#' # with restrictions:
#' check.gpt(file=file, latent=c("normal", "normal"),
#'              restrictions=list("g=.5","s1=s2=1", "m1_g=m1_d"))
#'              
#' # list of example model files:
#' list.files(path = paste0(path.package("gpt"),"/models"))
#' @export
check.gpt <- function(file, latent, restrictions=NULL){
  gpt <- new("gpt", file, latent, restrictions)
  gpt
}

