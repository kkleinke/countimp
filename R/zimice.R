# note: the function was imported from package mice 2.46.0 2017-10-23
zimice <-
function (data, m = 5, method = vector("character", length = ncol(data)),
          predictorMatrix = (1 - diag(1, ncol(data))), where = is.na(data),
          visitSequence = NULL, form = vector("character", length = ncol(data)),
          post = vector("character", length = ncol(data)), defaultMethod = c("pmm",
                                                                             "logreg", "polyreg", "polr"), maxit = 5, diagnostics = TRUE,
          printFlag = TRUE, seed = NA, imputationMethod = NULL, defaultImputationMethod = NULL,
          data.init = NULL, ...)
{
  call <- match.call()
  if (!is.na(seed))
    set.seed(seed)
  if (!(is.matrix(data) || is.data.frame(data)))
    stop("Data should be a matrix or data frame")
  nvar <- ncol(data)
  if (nvar < 2)
    stop("Data should contain at least two columns")
  data <- as.data.frame(data)
  nmis <- apply(is.na(data), 2, sum)
  varnames <- colnames(data)
  if (!(is.matrix(where) || is.data.frame(where)))
    stop("Argument `where` not a matrix or data frame")
  if (!all(dim(data) == dim(where)))
    stop("Arguments `data` and `where` not of same size")
  nwhere <- apply(where, 2, sum)
  state <- list(it = 0, im = 0, co = 0, dep = "", meth = "",
                log = FALSE)
  loggedEvents <- data.frame(it = 0, im = 0, co = 0, dep = "",
                             meth = "", out = "")
  if (!is.null(imputationMethod))
    method <- imputationMethod
  if (!is.null(defaultImputationMethod))
    defaultMethod <- defaultImputationMethod

  ## A NAMED `method` vector is matched by name, not by position.
  ##
  ## This engine was derived from mice 2.46, which stored `method` positionally
  ## and overwrote the user's names with the data's column names further down
  ## (names(method) <- varnames). The consequence was silent and severe: with
  ## data in the order (a, b), method = c(b = "poisson", a = "norm") applied
  ## "poisson" to `a` and "norm" to `b` -- the imputations came out under the
  ## wrong models, and nothing in the output said so, because the returned
  ## $method carried the data's names and therefore looked correct (B79).
  ## mice >= 3.0 matches by name; countimp now does too.
  ##
  ## An unnamed vector keeps positional semantics, which is what
  ## `method = "poisson"` (recycled) and `mice::mice()$method` rely on.
  if (!is.null(names(method)) && any(nzchar(names(method)))) {
    mn <- names(method)
    if (any(!nzchar(mn)))
      stop("`method` is partly named. Name every element or none, so that it ",
           "is unambiguous whether the vector is matched by name or by ",
           "position.", call. = FALSE)
    if (anyDuplicated(mn))
      stop("`method` names a variable more than once: ",
           paste(unique(mn[duplicated(mn)]), collapse = ", "), ".",
           call. = FALSE)
    unknown <- setdiff(mn, varnames)
    if (length(unknown))
      stop("`method` names variable(s) not in the data: ",
           paste(unknown, collapse = ", "),
           ".\n  Data columns: ", paste(varnames, collapse = ", "),
           call. = FALSE)
    ## Variables the user did not mention are marked NA, not "". The two mean
    ## different things and check.method() has to tell them apart:
    ##
    ##   ""  the user asked for this variable NOT to be imputed
    ##   NA  the user said nothing; use the default for the variable's type
    ##
    ## Conflating them lost data silently. check.method() only fills in defaults
    ## when *every* entry is "" (the mice 2.46 rule), so naming one variable left
    ## all the others unimputed -- and the named variable then came out
    ## incomplete too, because its own predictors still had missing values
    ## (B80: 17 of 60 cells in the named variable stayed NA).
    full <- setNames(rep(NA_character_, length(varnames)), varnames)
    full[mn] <- unname(method)
    method <- full
  }
  setup <- list(visitSequence = visitSequence, method = method,
                defaultMethod = defaultMethod, predictorMatrix = predictorMatrix,
                form = form, post = post, nvar = nvar, nmis = nmis, nwhere = nwhere,
                varnames = varnames)
  setup <- check.visitSequence(setup, where)
  setup <- check.method(setup, data)
  setup <- check.predictorMatrix(setup)
  setup <- check.data(setup, data, ...)
  .countimp_check_levels(data, setup$method, setup$predictorMatrix)
  method <- setup$method
  predictorMatrix <- setup$predictorMatrix
  visitSequence <- setup$visitSequence
  post <- setup$post
  p <- padModel(data, method, predictorMatrix, visitSequence,
                form, post, nvar)
  r <- (!is.na(p$data))
  imp <- vector("list", ncol(p$data))
  if (m > 0) {
    for (j in visitSequence) {
      y <- data[, j]
      ry <- r[, j]
      wy <- where[, j]
      imp[[j]] <- as.data.frame(matrix(NA, nrow = sum(wy),
                                       ncol = m))
      dimnames(imp[[j]]) <- list(row.names(data)[wy], 1:m)
      if (method[j] != "") {
        for (i in seq_len(m)) {
          if (nmis[j] < nrow(data)) {
            if (is.null(data.init)) {
              imp[[j]][, i] <- .countimp_impute_sample(y, ry, wy = wy)
            }
            else {
              imp[[j]][, i] <- data.init[wy, j]
            }
          }
          else imp[[j]][, i] <- rnorm(nrow(data))
        }
      }
    }
  }
  from <- 1
  to <- from + maxit - 1
  q <- zisampler(p, data, where, m, imp, r, visitSequence, c(from,
                                                           to), printFlag, ...)
  for (j in p$visitSequence) {
    p$data[!r[, j], j] <- NA
  }
  imp <- q$imp[seq_len(nvar)]
  names(imp) <- varnames
  names(method) <- varnames
  names(form) <- varnames
  names(post) <- varnames
  ## The engine indexes visitSequence numerically throughout (it is a copy of
  ## mice 2.46, where it was a column-number vector). mice >= 3.0 stores
  ## variable NAMES here, and mice.mids() relies on that. Convert at the
  ## boundary -- the engine keeps its integer version internally, the returned
  ## object carries the modern one.
  visitSequence <- varnames[visitSequence]
  if (!state$log)
    loggedEvents <- NULL
  if (state$log) {
    row.names(loggedEvents) <- seq_len(nrow(loggedEvents))
    ## mice 2.46 logged a `co` (column number) field that mice >= 3.0 dropped.
    ## mice.mids() rbind()s its own log onto this one and fails on the extra
    ## column, so drop it at the boundary. The information is redundant --
    ## `dep` names the same variable.
    loggedEvents$co <- NULL
  }
  ## `p$data` is padModel()'s working copy: it appends dummy columns and, for
  ## every factor used as a predictor, overwrites the contrast attribute with
  ## contr.treatment. Returning that copy as `$data` changed what the user got
  ## back: for an ORDERED factor R's default is contr.poly, so a model fitted
  ## on the completed data reported treatment contrasts (o2, o3) instead of
  ## polynomial ones (o.L, o.Q) -- different coefficients, silently, for data
  ## the user had passed in unmodified. The imputations themselves are
  ## unaffected (the internal coding is correct and deliberate); only the
  ## returned frame must carry the user's own attributes again.
  data_out <- as.data.frame(p$data[, seq_len(nvar)])
  for (v in intersect(names(data_out), names(data)))
    if (is.factor(data_out[[v]]) && is.factor(data[[v]]))
      attributes(data_out[[v]]) <- attributes(data[[v]])
  midsobj <- list(call = call, data = data_out,
                  where = where, m = m, nmis = nmis, imp = imp,
                  method = method, predictorMatrix = predictorMatrix, visitSequence = visitSequence,
                  form = form, post = post, seed = seed, iteration = q$iteration,
                  lastSeedValue = .Random.seed, chainMean = q$chainMean,
                  chainVar = q$chainVar, loggedEvents = loggedEvents, pad = p)

  ## Fields added by mice >= 3.0 -----------------------------------------------
  ## This engine is a copy of mice 2.46 (2017), so the object it built lacked
  ## the slots that mice 3.x added. Most downstream functions tolerate that
  ## (complete(), with(), pool(), plot() all worked), but mice.mids() reads
  ## `calltype` unconditionally and failed with an uninformative error about
  ## an argument of length zero. Filling the slots keeps countimp's objects
  ## usable with current mice without adopting mice's engine.
  varnames_all <- colnames(midsobj$data)
  midsobj$blocks   <- setNames(as.list(varnames_all), varnames_all)
  ## calltype "pred" tells mice that the models are defined by the
  ## predictorMatrix, not by formulas -- which is what this engine does.
  ## ("type" is not a value mice knows; it made mice.mids() look for a
  ## formula that is not there.)
  midsobj$calltype <- setNames(rep("pred", length(varnames_all)), varnames_all)
  ## mice keeps a formula per variable alongside the predictorMatrix and
  ## expects a real formula object in each slot, so derive them from the
  ## predictorMatrix rows rather than leaving them empty.
  midsobj$formulas <- setNames(lapply(varnames_all, function(v) {
    rhs <- varnames_all[predictorMatrix[v, ] != 0]
    stats::reformulate(if (length(rhs)) rhs else "1", response = v)
  }), varnames_all)
  midsobj$blots    <- setNames(vector("list", length(varnames_all)), varnames_all)
  midsobj$ignore   <- rep(FALSE, nrow(midsobj$data))
  midsobj$version  <- utils::packageVersion("countimp")
  midsobj$date     <- Sys.Date()

  if (!diagnostics)
    midsobj$pad <- NULL
  ## Own class first, mice's "mids" second (B55). Inheriting from "mids" keeps
  ## every mice function that tests inherits(x, "mids") working -- is.mids(),
  ## mice::complete(), mice::pool() -- while countimp's own S3 methods attach to
  ## countimp_mids and therefore never displace mice's methods for mice's own
  ## objects. Registering with.mids/print.mira/summary.mira directly, as earlier
  ## versions did, changed the behaviour of summary() on mice objects as a side
  ## effect of loading countimp.
  oldClass(midsobj) <- c("countimp_mids", "mids")
  return(midsobj)
}
