## ===========================================================================
## Repeated analysis of a mids object (since 3.0.0)
##
## countimp returns an object of class "mids" but shipped no with() method for
## it, relying on mice::with.mids. Without mice loaded, `with(imp, lm(y ~ x))`
## therefore dispatched to with.default(), which evaluates the expression in
## `imp` as a plain list -- so the analysis either failed with "object 'y' not
## found" or, when identically named variables happened to exist in the calling
## environment, ran SILENTLY on those instead of on the imputed data. In a
## simulation script, where the generating variables are usually still in the
## workspace, the second case is the normal one: the model was fitted to the
## complete cases only (n = 210 of 300 in the test case), the imputations were
## never used, and nothing indicated it.
##
## Providing the method here closes that hole and makes the analysis path
## independent of mice.
## ===========================================================================

#' Repeated Analysis of Multiply Imputed Data
#'
#' Evaluates an expression once per imputed data set.
#'
#' @param data An object of class \code{countimp_mids}, as returned by
#'   \code{\link{countimp}}. That class inherits from mice's \code{mids}, so
#'   mice functions that accept a \code{mids} object accept it too.
#' @param expr An expression to evaluate for each completed data set, e.g.
#'   \code{lm(y ~ x)}. Variables are looked up in the completed data first,
#'   then in the environment of the call.
#' @param ... Ignored, for compatibility.
#' @return An object of class \code{c("countimp_mira", "mira")} with components
#'   \code{call}, \code{call1}, \code{nmis} and \code{analyses} (a list of
#'   length \code{m}). The inherited \code{mira} class keeps the object usable
#'   with \code{mice::pool()}. Pass it to \code{\link{miinference}} to pool the
#'   results.
#' @details Compatible with \code{mice::with.mids} for the cases countimp
#'   supports, so existing analysis code keeps working. Unlike a plain
#'   \code{with()} on a list, the expression is evaluated against each of the
#'   \code{m} completed data sets, which is the entire point of multiple
#'   imputation.
#'
#'   The method is registered for \code{countimp_mids}, not for \code{mids}, so
#'   loading countimp leaves the behaviour of \code{with()}, \code{print()} and
#'   \code{summary()} on objects created by mice itself untouched. Objects
#'   created by \code{\link{countimp}} are dispatched here because
#'   \code{countimp_mids} precedes \code{mids} in their class vector.
#' @examples
#' data(crim4w)
#' imp <- countimp(crim4w, method = c(rep("", 5), "nb", "nb", "pmm", "pmm"),
#'                 m = 2, maxit = 1, printFlag = FALSE, seed = 1)
#' fit <- with(imp, lm(ACRIM ~ FEMALE + RE))
#' miinference(fit)
#' @author Kristian Kleinke

#' @export
with.countimp_mids <- function(data, expr, ...) {
  if (!inherits(data, "mids"))
    stop("`data` must be a mids object.", call. = FALSE)
  m    <- as.integer(data$m)
  call <- match.call()
  ## The expression is captured unevaluated and evaluated once per completed
  ## data set, with the caller's environment as the enclosure so that objects
  ## from the calling scope (weights, offsets, helper functions) stay visible.
  expr <- substitute(expr)
  env  <- parent.frame()

  analyses <- vector("list", m)
  for (i in seq_len(m)) {
    di <- countimp_complete(data, i)
    analyses[[i]] <- eval(expr, envir = di, enclos = env)
  }

  out <- list(call = call, call1 = data$call, nmis = data$nmis,
              analyses = analyses)
  ## Own class first, mice's "mira" second (B55) -- see the note in zimice.R.
  ## mice::mira() sets "mira" alone; the extra "matrix" this function used to
  ## add made inherits(x, "matrix") report TRUE while is.matrix(x) was FALSE and
  ## dim(x) was NULL, so code that branches on inherits() took the matrix path
  ## and then failed on the missing dimensions.
  oldClass(out) <- c("countimp_mira", "mira")
  out
}

#' @export
print.countimp_mira <- function(x, ...) {
  cat("Repeated analysis of", length(x$analyses), "imputed data sets\n")
  cat("Call: "); print(x$call)
  cat("\nFirst analysis:\n")
  print(x$analyses[[1L]])
  cat("\nUse miinference() to pool the", length(x$analyses), "analyses.\n")
  invisible(x)
}

#' @export
summary.countimp_mira <- function(object, ...) {
  lapply(object$analyses, summary, ...)
}
