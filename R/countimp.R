#' \code{countimp} Wrapper function
#' 
#' Depending on the entries of \code{mice}'s \code{method} argument, this function either calls \code{mice} directly or calls our variant to impute flat file zero-inflated count data.
#' 
#' @param ... \code{mice()} arguments see \code{help("mice")}
#' @author Kristian Kleinke
#' @examples 
#' data(crim4w)
#' ini <- countimp(crim4w, maxit=0)
#' meth <- ini$method
#' meth[6:7] <- "hp"
#' meth[8:9] <- "pmm"
#' pred <- ini$predictorMatrix
#' pred[,"id"] <- 0
#' pred["ACRIM",] <- c(0,1,3,2,0,3,3,2,1)
#' imp <- countimp( data = crim4w, method = meth, predictorMatrix = pred)
#' imp$method
#' @export
countimp <- function(...)
{
  arguments <- list(...)
  if (is.null(arguments$method)){do.call("mice",arguments)}else{
    if ( any(arguments$method %in% c("hnb", "hnb.boot", "hp", "hp.boot", "zinb", "zinb.boot", "zip", "zip.boot" ))  ){do.call("zimice",arguments)}else{do.call("mice",arguments)}}
}