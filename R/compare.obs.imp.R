#' Get summary statistics for the observed and imputed parts of the data
#' @param orig data frame of original (i.e. incomplete) data
#' @param imp list of length \code{m} of imputed data sets
#' @author Kristian Kleinke
#' @import mice
#' @importFrom stats aggregate
#' @export
compare.obs.imp <- function(orig,imp){
  nam <- names(imp$nmis[imp$nmis>0])
  output <- vector(length(nam), mode="list")
  names(output)<-nam
  for (i in 1:length(nam)){
    data <- mice::complete(imp,"long")
    ry <- is.na(orig[,nam[i]])
    ry[ry==TRUE] <- "imputed"
    ry[ry==FALSE] <- "observed"
    data <- cbind(data,ry)
    colnames(data)[length(colnames(data))] <- paste0("R.",nam[i])
    output[[i]] <- aggregate(data[,nam[i]], list(data[,paste0("R.",nam[i])]),summary)}
   return(output)
}
