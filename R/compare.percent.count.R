#' Compare percentages of observed and imputed counts
#' @param orig data frame of original (i.e. incomplete) data
#' @param imp list of length \code{m} of imputed data sets
#' @param counts vector giving the name(s) of the count variables for which observed and imputed counts shall be compared
#' @author Kristian Kleinke
#' @export
compare.percent.count <- function(orig,imp,counts){
  nam=names(orig[colSums(is.na(orig))>0])
  nam = nam[nam %in% counts]
  output <- vector(length(nam), mode="list")
  names(output)<-nam
  if (!length(nam))
    stop("None of `counts` is an incomplete variable in `orig`.", call. = FALSE)
  ## Stacked once, not once per variable (see compare.obs.imp).
  data0 <- .countimp_complete(imp, "long")
  for (i in 1:length(nam)){
    data <- data0
    ry= is.na(orig[,nam[i]])
    ry[ry==TRUE]<-"imputed"
    ry[ry==FALSE]<-"observed"
    data <- cbind(data,ry)
    colnames(data)[length(colnames(data))]<-paste0("R.",nam[i])
    tab <- table(data[,nam[i]],data[,paste0("R.",nam[i])])


    output[[i]]<- rbind("imputed"=round(tab[,1]/sum(tab[,1])*100,2),
                        "observed"=round(tab[,2]/sum(tab[,2])*100,2))

  }
  return(output)
}
