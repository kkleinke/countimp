#' @importFrom stats var
zisampler <-
function (p, data, where, m, imp, r, visitSequence, fromto, printFlag,
         ...)
{
  from <- fromto[1]
  to <- fromto[2]
  maxit <- to - from + 1
  chainVar <- chainMean <- NULL
  if (maxit > 0)
    chainVar <- chainMean <- array(0, dim = c(length(visitSequence),
                                              maxit, m), dimnames = list(dimnames(data)[[2]][visitSequence],
                                                                         seq_len(maxit), paste("Chain", seq_len(m))))
  if (maxit < 1)
    iteration <- 0
  if (maxit >= 1) {
    if (printFlag)
      cat("\n iter imp variable")
    for (k in from:to) {
      iteration <- k
      for (i in seq_len(m)) {
        if (printFlag)
          cat("\n ", iteration, " ", i)
        for (j in visitSequence) {
          wy <- where[, j]
          ry <- r[, j]
          p$data[(!ry) & wy, j] <- imp[[j]][(!ry)[wy],
                                            i]
        }
        for (j in setdiff(p$visitSequence, visitSequence)) {
          cat.columns <- p$data[, p$categories[j, 4]]
          mm <- model.matrix(~cat.columns - 1, model.frame(~cat.columns,
                                                           na.action = na.pass))[, -1]
          p$data[, (j:(j + p$categories[p$categories[j,
                                                     4], 2] - 1))] <- mm
        }
        for (j in p$visitSequence) {
          theMethod <- p$method[j]
          vname <- dimnames(p$data)[[2]][j]
          oldstate <- get("state", pos = parent.frame())
          newstate <- list(it = k, im = i, co = j, dep = vname,
                           meth = theMethod, log = oldstate$log)
          assign("state", newstate, pos = parent.frame(),
                 inherits = TRUE)
          if (printFlag && theMethod != "dummy")
            cat(" ", vname)
          empt <- theMethod == ""
          elem <- !empt && !is.passive(theMethod) &&
            theMethod != "dummy"
          flat <- elem && substring(theMethod, 1, 2) !=
            "2l"
          pass <- !empt && is.passive(theMethod)
          dumm <- theMethod == "dummy"
          if (elem) {
            if (flat) {
              predictors <- p$predictorMatrix[j, ] !=
                0
            }
            else {
              predictors <- p$predictorMatrix[j, ] !=
                0
            }
            if (!is.null(p$form) && nchar(p$form[j]) >
                0) {
              myform <- paste(p$form[j], "0", sep = "+")
              x <- model.matrix(formula(myform), p$data)
            }
            else {
              x <- p$data[, predictors, drop = FALSE]
            }
            y <- p$data[, j]
            ry <- complete.cases(x, y) & r[, j]
            wy <- complete.cases(x) & where[, j]
            cc <- wy[where[, j]]
            type <- p$predictorMatrix[j, predictors]
            if (k == 1)
              check.df(x, y, ry)
            keep <- remove.lindep(x, y, ry, ...)
            x <- x[, keep, drop = FALSE]
            type <- type[keep]
            f <- paste("mice.impute", theMethod, sep = ".")
            imputes <- p$data[wy, j]
            imputes[!cc] <- NA
            imputes[cc] <- do.call(f, args = list(y,
                                                  ry, x, wy = wy, type = type, ...))
            imp[[j]][, i] <- imputes
            p$data[(!r[, j]) & where[, j], j] <- imp[[j]][(!r[,
                                                              j])[where[, j]], i]
          }
          if (pass) {
            wy <- where[, j]
            imp[[j]][, i] <- model.frame(as.formula(theMethod),
                                         p$data[wy, ], na.action = na.pass)
            p$data[(!ry) & wy, j] <- imp[[j]][(!ry)[wy],
                                              i]
          }
          if (dumm) {
            cat.columns <- p$data[, p$categories[j, 4]]
            mm <- model.matrix(~cat.columns - 1, model.frame(~cat.columns,
                                                             na.action = na.pass))[, -1]
            p$data[, (j:(j + p$categories[p$categories[j,
                                                       4], 2] - 1))] <- mm
            remove("cat.columns")
          }
          cmd <- p$post[j]
          if (cmd != "") {
            eval(parse(text = cmd))
            p$data[where[, j], j] <- imp[[j]][, i]
          }
        }
      }
      k2 <- k - from + 1
      if (length(visitSequence) > 0) {
        for (j in seq_along(visitSequence)) {
          jj <- visitSequence[j]
          if (!is.factor(data[, jj])) {
            chainVar[j, k2, ] <- apply(imp[[jj]], 2,
                                       var, na.rm = TRUE)
            chainMean[j, k2, ] <- colMeans(as.matrix(imp[[jj]]),
                                           na.rm = TRUE)
          }
          if (is.factor(data[, jj])) {
            for (mm in seq_len(m)) {
              nc <- as.integer(factor(imp[[jj]][, mm],
                                      levels = levels(data[, jj])))
              chainVar[j, k2, mm] <- var(nc, na.rm = TRUE)
              chainMean[j, k2, mm] <- mean(nc, na.rm = TRUE)
            }
          }
        }
      }
    }
    if (printFlag)
      cat("\n")
  }
  return(list(iteration = maxit, imp = imp, chainMean = chainMean,
              chainVar = chainVar))
}

is.passive <-
function (string) 
{
  return("~" == substring(string, 1, 1))
}

check.df <-
function (x, y, ry) 
{
  df <- sum(ry) - ncol(x) - 1
  mess <- paste("df set to 1. # observed cases:", sum(ry), 
                " # predictors:", ncol(x) + 1)
  if (df < 1) 
    updateLog(out = mess, frame = 3)
}

remove.lindep <-
function (x, y, ry, eps = 1e-04, maxcor = 0.99, allow.na = FALSE, 
          ...) 
{
  if (ncol(x) == 0) 
    return(NULL)
  if (eps <= 0) 
    stop("\n Argument 'eps' must be positive.")
  if (allow.na && sum(ry) == 0) {
    updateLog(out = "No observed outcomes, keep all predictors", 
              frame = 3)
    return(rep.int(TRUE, ncol(x)))
  }
  xobs <- x[ry, , drop = FALSE]
  yobs <- as.numeric(y[ry])
  keep <- unlist(apply(xobs, 2, var) > eps)
  keep[is.na(keep)] <- FALSE
  highcor <- suppressWarnings(unlist(apply(xobs, 2, cor, yobs) < 
                                       maxcor))
  keep <- keep & highcor
  if (all(!keep)) 
    updateLog(out = "All predictors are constant or have too high correlation.", 
              frame = 3)
  if (length(keep) == 1) 
    keep[1] <- TRUE
  k <- sum(keep)
  cx <- cor(xobs[, keep, drop = FALSE], use = "all.obs")
  eig <- eigen(cx, symmetric = TRUE)
  ncx <- cx
  while (eig$values[k]/eig$values[1] < eps) {
    j <- seq_len(k)[order(abs(eig$vectors[, k]), decreasing = TRUE)[1]]
    keep[keep][j] <- FALSE
    ncx <- cx[keep[keep], keep[keep], drop = FALSE]
    k <- k - 1
    eig <- eigen(ncx)
  }
  if (!all(keep)) {
    out <- paste(dimnames(x)[[2]][!keep], collapse = ", ")
    updateLog(out = out, frame = 3)
  }
  return(keep)
}

updateLog <-
function (out = NULL, meth = NULL, frame = 2) 
{
  s <- get("state", parent.frame(frame))
  r <- get("loggedEvents", parent.frame(frame))
  rec <- data.frame(it = s$it, im = s$im, co = s$co, dep = s$dep, 
                    meth = if (is.null(meth)) 
                      s$meth
                    else meth, out = if (is.null(out)) 
                      ""
                    else out)
  if (s$log) 
    rec <- rbind(r, rec)
  s$log <- TRUE
  assign("state", s, pos = parent.frame(frame), inherits = TRUE)
  assign("loggedEvents", rec, pos = parent.frame(frame), inherits = TRUE)
  return()
}

check.visitSequence <-
function (setup, where) 
{
  nwhere <- setup$nwhere
  nvar <- setup$nvar
  visitSequence <- setup$visitSequence
  if (is.null(visitSequence)) 
    visitSequence <- seq_len(ncol(where))[apply(where, 2, 
                                                any)]
  if (!is.numeric(visitSequence)) {
    code <- match.arg(visitSequence, c("roman", "arabic", 
                                       "monotone", "revmonotone"))
    visitSequence <- switch(code, roman = seq_len(nvar)[nwhere > 
                                                          0], arabic = rev(seq_len(nvar)[nwhere > 0]), monotone = order(nwhere)[(sum(nwhere == 
                                                                                                                                       0) + 1):length(nwhere)], revmonotone = rev(order(nwhere)[(sum(nwhere == 
                                                                                                                                                                                                       0) + 1):length(nwhere)]), seq_len(nvar)[nwhere > 
                                                                                                                                                                                                                                                 0])
  }
  flags <- nwhere == 0 & is.element(seq_len(nvar), visitSequence)
  if (any(flags)) 
    visitSequence <- visitSequence[!flags]
  visitSequence <- visitSequence[visitSequence <= nvar]
  visitSequence <- visitSequence[visitSequence >= 1]
  setup$visitSequence <- visitSequence
  return(setup)
}

check.method <-
function (setup, data) 
{
  method <- setup$method
  defaultMethod <- setup$defaultMethod
  visitSequence <- setup$visitSequence
  nwhere <- setup$nwhere
  nvar <- setup$nvar
  if (all(method == "")) {
    for (j in visitSequence) {
      y <- data[, j]
      if (is.numeric(y)) {
        method[j] <- defaultMethod[1]
      }
      else if (nlevels(y) == 2) {
        method[j] <- defaultMethod[2]
      }
      else if (is.ordered(y) && nlevels(y) > 2) {
        method[j] <- defaultMethod[4]
      }
      else if (nlevels(y) > 2) {
        method[j] <- defaultMethod[3]
      }
      else if (is.logical(y)) {
        method[j] <- defaultMethod[2]
      }
      else {
        method[j] <- defaultMethod[1]
      }
    }
  }
  if (length(method) == 1) {
    if (is.passive(method)) 
      stop("Cannot have a passive imputation method for every column.")
    method <- rep(method, nvar)
  }
  if (length(method) != nvar) {
    stop(paste0("The length of method (", length(method), 
                ") does not match the number of columns in the data (", 
                nvar, ")."))
  }
  active.check <- !is.passive(method) & nwhere > 0 & method != 
    ""
  passive.check <- is.passive(method) & nwhere > 0 & method != 
    ""
  check <- all(active.check) & any(passive.check)
  if (check) {
    fullNames <- rep.int("mice.impute.passive", length(method[passive.check]))
  }
  else {
    fullNames <- paste("mice.impute", method[active.check], 
                       sep = ".")
    if (length(method[active.check]) == 0) 
      fullNames <- character(0)
  }
  notFound <- !vapply(fullNames, exists, logical(1), mode = "function", 
                      inherits = TRUE)
  if (any(notFound)) {
    stop(paste("The following functions were not found:", 
               paste(fullNames[notFound], collapse = ", ")))
  }
  for (j in visitSequence) {
    y <- data[, j]
    vname <- dimnames(data)[[2]][j]
    mj <- method[j]
    mlist <- list(m1 = c("logreg", "logreg.boot", "polyreg", 
                         "lda", "polr"), m2 = c("norm", "norm.nob", "norm.predict", 
                                                "norm.boot", "mean", "2l.norm", "2l.pan", "2lonly.pan", 
                                                "quadratic", "ri"), m3 = c("norm", "norm.nob", "norm.predict", 
                                                                           "norm.boot", "mean", "2l.norm", "2l.pan", "2lonly.pan", 
                                                                           "quadratic", "logreg", "logreg.boot"))
    if (is.numeric(y) && (mj %in% mlist$m1)) {
      warning("Type mismatch for variable ", vname, "\nImputation method ", 
              mj, " is for categorical data.", "\nIf you want that, turn variable ", 
              vname, " into a factor,", "\nand store your data in a data frame.", 
              call. = FALSE)
    }
    else if (is.factor(y) && nlevels(y) == 2 && (mj %in% 
                                                 mlist$m2)) {
      warning("Type mismatch for variable ", vname, "\nImputation method ", 
              mj, " is not for factors.", call. = FALSE)
    }
    else if (is.factor(y) && nlevels(y) > 2 && (mj %in% mlist$m3)) {
      warning("Type mismatch for variable ", vname, "\nImputation method ", 
              mj, " is not for factors with three or more levels.", 
              call. = FALSE)
    }
  }
  setup$method <- method
  return(setup)
}

check.predictorMatrix <-
function (setup) 
{
  pred <- setup$predictorMatrix
  varnames <- setup$varnames
  nwhere <- setup$nwhere
  nvar <- setup$nvar
  vis <- setup$visitSequence
  post <- setup$post
  if (!is.matrix(pred)) 
    stop("Argument 'predictorMatrix' not a matrix.")
  if (nvar != nrow(pred) || nvar != ncol(pred)) 
    stop(paste("The predictorMatrix has", nrow(pred), "rows and", 
               ncol(pred), "columns. Both should be", nvar, "."))
  dimnames(pred) <- list(varnames, varnames)
  diag(pred) <- 0
  for (j in seq_len(nvar)) {
    if (nwhere[j] == 0 && any(pred[j, ] != 0)) 
      pred[j, ] <- 0
  }
  setup$predictorMatrix <- pred
  setup$visitSequence <- vis
  setup$post <- post
  return(setup)
}

check.data <-
function (setup, data, allow.na = FALSE, remove_collinear = TRUE, 
          ...) 
{
  pred <- setup$predictorMatrix
  nvar <- setup$nvar
  varnames <- setup$varnames
  meth <- setup$method
  vis <- setup$visitSequence
  post <- setup$post
  isclassvar <- apply(pred == -2, 2, any)
  for (j in seq_len(nvar)) {
    if (isclassvar[j] && is.factor(data[, j])) 
      stop("Class variable (column ", j, ") cannot be factor. Convert to numeric by as.integer()")
  }
  for (j in seq_len(nvar)) {
    if (!is.passive(meth[j])) {
      d.j <- data[, j]
      v <- if (is.character(d.j)) 
        NA
      else var(as.numeric(d.j), na.rm = TRUE)
      constant <- if (allow.na) {
        if (is.na(v)) 
          FALSE
        else v < 1000 * .Machine$double.eps
      }
      else {
        is.na(v) || v < 1000 * .Machine$double.eps
      }
      didlog <- FALSE
      if (constant && any(pred[, j] != 0)) {
        out <- varnames[j]
        pred[, j] <- 0
        updateLog(out = out, meth = "constant")
        didlog <- TRUE
      }
      if (constant && meth[j] != "") {
        out <- varnames[j]
        pred[j, ] <- 0
        if (!didlog) 
          updateLog(out = out, meth = "constant")
        meth[j] <- ""
        vis <- vis[vis != j]
        post[j] <- ""
      }
    }
  }
  ispredictor <- apply(pred != 0, 2, any)
  if (any(ispredictor)) {
    droplist <- find.collinear(data[, ispredictor, drop = FALSE], 
                               ...)
  }
  else {
    droplist <- NULL
  }
  if (length(droplist) > 0 && remove_collinear) {
    for (k in seq_along(droplist)) {
      j <- which(varnames %in% droplist[k])
      didlog <- FALSE
      if (any(pred[, j] != 0)) {
        out <- varnames[j]
        pred[, j] <- 0
        updateLog(out = out, meth = "collinear")
        didlog <- TRUE
      }
      if (meth[j] != "") {
        out <- varnames[j]
        pred[j, ] <- 0
        if (!didlog) 
          updateLog(out = out, meth = "collinear")
        meth[j] <- ""
        vis <- vis[vis != j]
        post[j] <- ""
      }
    }
  }
  setup$predictorMatrix <- pred
  setup$visitSequence <- vis
  setup$post <- post
  setup$meth <- meth
  return(setup)
}

padModel <-
function (data, method, predictorMatrix, visitSequence, form, 
          post, nvar) 
{
  categories <- data.frame(is.factor = factor(rep.int(FALSE, 
                                                      nvar), levels = c("TRUE", "FALSE")), n.dummy = rep.int(0, 
                                                                                                             nvar), is.dummy = factor(rep.int(FALSE, nvar), levels = c("TRUE", 
                                                                                                                                                                       "FALSE")), father = rep.int(0, nvar))
  data <- data
  for (j in seq_len(nvar)) {
    if (is.factor(data[, j]) && any(predictorMatrix[, j] != 
                                    0)) {
      categories[j, 1] <- TRUE
      data[, j] <- C(data[, j], contr.treatment)
      n.dummy <- length(levels(data[, j])) - 1
      categories[j, 2] <- n.dummy
      predictorMatrix <- rbind(predictorMatrix, matrix(0, 
                                                       ncol = ncol(predictorMatrix), nrow = n.dummy))
      predictorMatrix <- cbind(predictorMatrix, matrix(rep(predictorMatrix[, 
                                                                           j], times = n.dummy), ncol = n.dummy))
      predictorMatrix[seq_len(nvar), j] <- rep.int(0, times = nvar)
      form <- c(form, rep.int("", n.dummy))
      if (any(visitSequence == j)) {
        idx <- (ncol(predictorMatrix) - n.dummy + 1):ncol(predictorMatrix)
        predictorMatrix[idx, j] <- rep.int(1, times = n.dummy)
        newcol <- ncol(predictorMatrix) - n.dummy + 1
        nloops <- sum(visitSequence == j)
        for (ii in seq_len(nloops)) {
          idx2 <- seq_along(visitSequence)[visitSequence == 
                                             j][ii]
          visitSequence <- append(visitSequence, newcol, 
                                  idx2)
        }
      }
      data <- cbind(data, matrix(0, ncol = n.dummy, nrow = nrow(data)))
      idx <- (ncol(predictorMatrix) - n.dummy + 1):ncol(predictorMatrix)
      data[is.na(data[, j]), idx] <- NA
      cat.column <- data[!is.na(data[, j]), j]
      data[!is.na(data[, j]), idx] <- model.matrix(~cat.column - 
                                                     1)[, -1]
      names(data)[idx] <- paste(attr(data, "names")[j], 
                                seq_len(n.dummy), sep = ".")
      method <- c(method, rep.int("dummy", n.dummy))
      post <- c(post, rep.int("", n.dummy))
      categories <- rbind(categories, data.frame(is.factor = rep.int(FALSE, 
                                                                     n.dummy), n.dummy = rep.int(0, n.dummy), is.dummy = rep.int(TRUE, 
                                                                                                                                 n.dummy), father = rep.int(j, n.dummy)))
    }
  }
  varnames <- dimnames(data)[[2]]
  dimnames(predictorMatrix) <- list(varnames, varnames)
  names(method) <- varnames
  names(form) <- varnames
  names(post) <- varnames
  names(visitSequence) <- varnames[visitSequence]
  dimnames(categories)[[1]] <- dimnames(data)[[2]]
  if (anyDuplicated(names(data))) 
    stop("Column names of padded data not unique")
  return(list(data = as.data.frame(data), predictorMatrix = predictorMatrix, 
              method = method, visitSequence = visitSequence, form = form, 
              post = post, categories = categories))
}

find.collinear <-
function (x, threshold = 0.999, ...) 
{
  nvar <- ncol(x)
  x <- data.matrix(x)
  r <- !is.na(x)
  nr <- apply(r, 2, sum, na.rm = TRUE)
  ord <- order(nr, decreasing = TRUE)
  xo <- x[, ord, drop = FALSE]
  varnames <- dimnames(xo)[[2]]
  z <- suppressWarnings(cor(xo, use = "pairwise.complete.obs"))
  hit <- outer(seq_len(nvar), seq_len(nvar), "<") & (abs(z) >= 
                                                       threshold)
  out <- apply(hit, 2, any, na.rm = TRUE)
  return(varnames[out])
}