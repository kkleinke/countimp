## ===========================================================================
## Extracting completed data sets from a mids object (since 3.0.0)
##
## This was the last hard call into mice inside countimp: compare.obs.imp()
## and compare.percent.count() used mice::complete(). A mids object stores the
## incomplete data plus, per variable, an m-column data frame of imputations
## with the row names giving the target rows -- filling that in is a few lines,
## and doing it here removes countimp's dependency on the layout of mice's
## accessor rather than on the object it returns.
##
## The reader is deliberately tolerant: it works on countimp's own mids
## objects and on mice's, and it does not assume that `where` is present
## (mice added it in 3.0; objects from 2.x do not carry it).
## ===========================================================================

#' Extract Completed Data from an Imputation Object
#'
#' Fills the missing values of the original data with the imputations stored in
#' a \code{mids} object.
#'
#' @param data An object of class \code{mids}, as returned by
#'   \code{\link{countimp}}.
#' @param action Which completed data to return:
#'   \describe{
#'     \item{an integer \code{1..m}}{that single completed data set.}
#'     \item{\code{"long"}}{all \code{m} data sets stacked, with leading
#'       columns \code{.imp} (imputation number) and \code{.id} (row name).}
#'     \item{\code{"all"} or \code{"list"}}{a list of length \code{m}.}
#'     \item{\code{"broad"}}{all \code{m} data sets side by side, columns
#'       suffixed with \code{.1 ... .m}.}
#'   }
#' @param include Logical; prepend the original incomplete data as imputation
#'   \code{0}. Only used for \code{"long"}, \code{"all"} and \code{"broad"}.
#' @param ... Ignored, for compatibility.
#' @return A data frame, or a list of data frames for \code{action = "all"}.
#' @details Written to be interchangeable with \code{mice::complete()} for the
#'   actions listed above, so that code using either works unchanged. countimp
#'   uses this function internally, which is what makes its post-imputation
#'   tools independent of \pkg{mice} being installed.
#' @examples
#' data(crim4w)
#' imp <- countimp(crim4w, method = c(rep("", 5), "nb", "nb", "pmm", "pmm"),
#'                 m = 2, maxit = 1, printFlag = FALSE, seed = 1)
#' dim(countimp_complete(imp, 1))
#' dim(countimp_complete(imp, "long"))
#' @author Kristian Kleinke
#' @export
countimp_complete <- function(data, action = 1L, include = FALSE, ...) {
  if (!inherits(data, "mids"))
    stop("`data` must be a mids object, as returned by countimp().",
         call. = FALSE)
  m    <- as.integer(data$m)
  orig <- data$data

  ## Which cells were imputed. `where` is authoritative when present; older
  ## objects only record the missing-data pattern.
  wh <- data$where
  if (is.null(wh)) wh <- is.na(orig)

  single <- function(i) {
    out <- orig
    for (v in names(data$imp)) {
      im <- data$imp[[v]]
      if (is.null(im) || !nrow(im) || !(v %in% names(out))) next
      rows <- rownames(im)
      idx  <- if (!is.null(rows) && all(rows %in% rownames(orig)))
                match(rows, rownames(orig)) else which(wh[, v])
      if (!length(idx)) next
      val <- im[[i]]
      ## Keep the column's type: a factor column must stay a factor with its
      ## original levels even if the imputations arrive as characters.
      if (is.factor(out[[v]]) && !is.factor(val))
        val <- factor(as.character(val), levels = levels(out[[v]]))
      out[idx, v] <- val
    }
    out
  }

  if (is.numeric(action) && length(action) == 1L) {
    i <- as.integer(action)
    if (i == 0L) return(orig)
    if (i < 1L || i > m)
      stop("`action` must be between 1 and m = ", m, ".", call. = FALSE)
    return(single(i))
  }

  action <- match.arg(as.character(action),
                      c("long", "all", "list", "broad", "repeated"))
  idx <- if (isTRUE(include)) 0L:m else seq_len(m)
  sets <- lapply(idx, function(i) if (i == 0L) orig else single(i))

  if (action %in% c("all", "list")) {
    names(sets) <- as.character(idx)
    ## mice::complete(action = "all") returns the list with class "mild". mice
    ## itself defines no methods for it, but user code and downstream packages
    ## do test inherits(x, "mild"), so a drop-in replacement has to carry it --
    ## B10 caught the difference. "list" stays second so that every list
    ## operation keeps working unchanged.
    return(structure(sets, class = c("mild", "list")))
  }

  if (action %in% c("broad", "repeated")) {
    out <- do.call(cbind, sets)
    names(out) <- paste0(rep(names(orig), length(idx)), ".",
                         rep(idx, each = ncol(orig)))
    return(out)
  }

  ## long. Column order and the type of `.id` follow mice::complete(): the
  ## data columns come first, `.imp`/`.id` last, and `.id` is integer when the
  ## row names are the default integer sequence. Code that subsets by position
  ## or feeds `.id` to something numeric therefore behaves the same.
  out <- do.call(rbind, sets)
  rn  <- rownames(orig)
  id  <- if (identical(rn, as.character(seq_len(nrow(orig))))) seq_len(nrow(orig)) else rn
  cbind(out,
        .imp = rep(idx, each = nrow(orig)),
        .id  = rep(id, length(idx)),
        row.names = NULL)
}

## Internal accessor: prefer countimp's own reader, fall back to mice's for
## objects that countimp cannot read (should not occur, but a mids object from
## a future mice version might carry a layout this function does not expect).
.countimp_complete <- function(imp, action = "long") {
  out <- try(countimp_complete(imp, action), silent = TRUE)
  if (!inherits(out, "try-error")) return(out)
  .countimp_need_mice("reading this mids object")
  mice::complete(imp, action)
}
