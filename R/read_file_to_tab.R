
###############################
# file: path to table

#' @importFrom utils read.delim
read.file.to.tab <- function(file){
  
  if (grepl("\n", file)){
    file <- textConnection(file)
  } 
  tab <- read.delim(file, sep=";", 
                    header = FALSE, strip.white = TRUE, 
                    stringsAsFactors = FALSE, blank.lines.skip = TRUE, 
                    comment.char="#")
  
  
  # number of continuous variables:
  n.cont <- ncol(tab) - 3
  if (n.cont<1)
    stop ("Specification of continuous distributions missing in model file.")
  
  # remove spaces
  for (i in 1:ncol(tab)){
    tab[,i] <- gsub(" ","",tab[,i], fixed = TRUE)
  }
  if (any(tab == ""))
    stop ("Model file not completely specified (missing/empty entries)!")
  
  colnames(tab) <- c("tree", "cat", "mpt", paste0("c",1:n.cont))
  
  tab
}
