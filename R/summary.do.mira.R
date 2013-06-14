summary.do.mira <- function(object, ...)
{
  cat("\n Pooled Fixed Effects Coefficients:\n")
  printCoefmat(as.matrix(object$MIINFERENCE[,1:5]), digits=3, P.value=TRUE, has.Pvalue=TRUE, signif.stars=TRUE,signif.legend=TRUE)
  print(object$MIINFERENCE[,6:9])
  cat("\n Pooled Random Effects SD(s):\n")
  print(object$PR.SD)
  cat("\n Pooled Random Effects Correlation(s):\n")
  print(object$PR.COR)
}