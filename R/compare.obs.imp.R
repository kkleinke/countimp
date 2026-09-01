#' Get summary statistics for the observed and imputed parts of the data
#' @param orig data frame of original (i.e. incomplete) data
#' @param imp list of length \code{m} of imputed data sets
#' @author Kristian Kleinke
#' @importFrom stats aggregate
#' @export
compare.obs.imp <- function(orig,imp){
  nam <- names(imp$nmis[imp$nmis>0])
  if (!length(nam))
    stop("No incomplete variables in this imputation object.", call. = FALSE)
  output <- vector(length(nam), mode="list")
  names(output)<-nam
  ## The stacked data are the same for every variable; building them inside the
  ## loop repeated the work once per incomplete variable. Reading them through
  ## countimp's own accessor also removes the dependency on mice being
  ## installed for a purely descriptive summary.
  data0 <- .countimp_complete(imp, "long")
  for (i in 1:length(nam)){
    data <- data0
    ry <- is.na(orig[,nam[i]])
    ry[ry==TRUE] <- "imputed"
    ry[ry==FALSE] <- "observed"
    data <- cbind(data,ry)
    colnames(data)[length(colnames(data))] <- paste0("R.",nam[i])
    output[[i]] <- aggregate(data[,nam[i]], list(data[,paste0("R.",nam[i])]),summary)}
   return(output)
}
