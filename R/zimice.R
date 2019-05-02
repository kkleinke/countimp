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
  setup <- list(visitSequence = visitSequence, method = method,
                defaultMethod = defaultMethod, predictorMatrix = predictorMatrix,
                form = form, post = post, nvar = nvar, nmis = nmis, nwhere = nwhere,
                varnames = varnames)
  setup <- check.visitSequence(setup, where)
  setup <- check.method(setup, data)
  setup <- check.predictorMatrix(setup)
  setup <- check.data(setup, data, ...)
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
              imp[[j]][, i] <- mice::mice.impute.sample(y,
                                                  ry, wy = wy, ...)
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
  names(visitSequence) <- varnames[visitSequence]
  if (!state$log)
    loggedEvents <- NULL
  if (state$log)
    row.names(loggedEvents) <- seq_len(nrow(loggedEvents))
  midsobj <- list(call = call, data = as.data.frame(p$data[,
                                                           seq_len(nvar)]), where = where, m = m, nmis = nmis, imp = imp,
                  method = method, predictorMatrix = predictorMatrix, visitSequence = visitSequence,
                  form = form, post = post, seed = seed, iteration = q$iteration,
                  lastSeedValue = .Random.seed, chainMean = q$chainMean,
                  chainVar = q$chainVar, loggedEvents = loggedEvents, pad = p)
  if (!diagnostics)
    midsobj$pad <- NULL
  oldClass(midsobj) <- "mids"
  return(midsobj)
}
