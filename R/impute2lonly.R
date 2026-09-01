## impute2lonly.R -- variables that live at a higher level than the rows
##
## Why this file exists. In a two- or three-level data set some variables do not
## vary within their cluster: class size, school type, a therapist's years of
## experience. In the long format they still occupy one cell per row, and the
## FCS engine imputes cell by cell -- so an incomplete class-level variable came
## back with a different value in every row of the same class.
##
## Measured before this was built (25 replications, 60 classes, 20% of the
## CLASSES missing their class variable, analyse/k28_ebenen_belege.R):
##
##   classes whose constancy was destroyed   12 of 60
##   effect of the class variable            0.486 against 0.507 on complete
##                                           data, true value 0.500
##
## The coefficient survives in that design, which is exactly what makes this
## dangerous: nothing in the output says that the completed data now contain
## classes whose class size differs between their own pupils. Any analysis that
## aggregates to the class level reads that as real variation.
##
## The remedy is the one mice uses for the same problem (2lonly.norm,
## 2lonly.pmm): impute at the level the variable lives at. Aggregate to one row
## per cluster, impute there, and hand the result back to every row of the
## cluster. countimp adds the count families to that list, because a
## cluster-level count -- beds on a ward, staff in a practice -- is the case
## this package exists for.


## Is `y` constant within each level of `g`, over the OBSERVED values?
##
## `ry` marks what is observed and must be passed wherever it is known. Inside
## an imputation run, is.na() is NOT a substitute: the engine fills every
## missing cell with a hot-deck starting value before the first draw, so a
## cluster-level variable arrives row-wise varying and NOT missing. Measured:
## without ry the first call refused its own data with "not constant within any
## of the grouping variable(s)". Outside a run (the setup check below) the data
## still carry NA, and ry = NULL is right there.
.countimp_constant_within <- function(y, g, ry = NULL) {
  if (!length(y)) return(TRUE)
  ok <- if (is.null(ry)) !is.na(y) else ry & !is.na(y)
  if (!any(ok)) return(TRUE)
  n <- tapply(y[ok], g[ok], function(z) length(unique(z)))
  all(n <= 1L, na.rm = TRUE)
}


## Which grouping variable does a 2lonly method work on?
##
## Same convention as every other two-level method: the column coded -2. With
## several (three-level data) the variable may live at any of them, so the
## OUTERMOST level it is constant within is taken -- a school-level variable is
## constant within classes too, and imputing it per class would produce
## different values for classes of the same school.
.countimp_2lonly_grp <- function(y, x, type, ry = NULL) {
  grp <- which(type == -2L)
  if (!length(grp))
    stop("countimp: a 2lonly method needs a grouping variable: code one ",
         "predictor as -2 in `type`.", call. = FALSE)
  nam <- colnames(x)
  kand <- nam[grp]
  konst <- vapply(kand, function(g) .countimp_constant_within(y, x[, g], ry),
                  logical(1))
  if (!any(konst))
    stop("countimp: the variable to impute is not constant within any of the ",
         "grouping\n  variable(s) ", paste(kand, collapse = ", "),
         ", so it does not live at a cluster level.\n",
         "  A 2lonly method is the wrong choice here; use the ordinary method ",
         "for its\n  measurement level.", call. = FALSE)
  ## outermost = fewest distinct values among those it is constant within
  moeg <- kand[konst]
  moeg[which.min(vapply(moeg, function(g) length(unique(x[, g])), integer(1)))]
}


## The engine. `basis` is the name of a countimp.impute.* method that does the
## work on the aggregated data.
##
## Aggregation. One row per cluster: the variable itself takes its (unique)
## observed value, the predictors take their cluster mean. Averaging the
## predictors is what mice's 2lonly.* do, and it is the only choice that does
## not depend on row order -- a level-2 predictor is already constant, so its
## mean IS its value, and a level-1 predictor enters as the cluster mean, which
## is the quantity a cluster-level model can use.
##
## The grouping columns themselves are dropped from the predictor matrix: their
## identifiers are labels, not covariates, and passing them on would let the
## imputation model regress on cluster numbering.
.countimp_2lonly <- function(y, ry, x, type, basis, wy = NULL, ...) {
  if (is.null(wy)) wy <- !ry
  x <- as.data.frame(x)
  g.name <- .countimp_2lonly_grp(y, x, type, ry)
  g <- as.character(x[[g.name]])

  ## Cluster level: one row per cluster, built from the OBSERVED cells only.
  ## The starting values the engine wrote into the missing cells are ordinary
  ## numbers, so !is.na(y) would take them for data and there would be nothing
  ## left to impute.
  cl <- unique(g)
  yc <- vapply(cl, function(k) {
    v <- y[g == k & ry & !is.na(y)]
    if (length(v)) as.numeric(v[1L]) else NA_real_
  }, numeric(1))

  grp.cols <- colnames(x)[type == -2L]
  pred.cols <- setdiff(colnames(x), grp.cols)
  xc <- if (length(pred.cols)) {
    m <- vapply(pred.cols, function(p)
      vapply(cl, function(k) mean(as.numeric(x[[p]])[g == k], na.rm = TRUE),
             numeric(1)), numeric(length(cl)))
    matrix(m, nrow = length(cl), dimnames = list(NULL, pred.cols))
  } else matrix(1, nrow = length(cl), ncol = 1L,
                dimnames = list(NULL, "konstante"))
  ## A predictor that is missing for a whole cluster becomes NaN in the mean.
  ## Left as is it would silently drop that cluster from the fit; replaced by
  ## the overall mean it is at least a defined, documented fallback.
  for (j in seq_len(ncol(xc))) {
    schlecht <- !is.finite(xc[, j])
    if (any(schlecht)) {
      xc[schlecht, j] <- mean(xc[!schlecht, j])
      .countimp_note_event("2lonly_predictor_all_missing",
        sprintf("predictor '%s' missing for %d whole cluster(s); cluster mean replaced by the overall mean",
                colnames(xc)[j], sum(schlecht)))
    }
  }

  ryc <- !is.na(yc)
  if (!any(ryc))
    stop("countimp: the variable is missing in every cluster, so there is ",
         "nothing to impute from.", call. = FALSE)
  if (all(ryc)) return(y[wy])          # nothing missing at cluster level

  f <- .countimp_find_method(basis)
  if (is.null(f))
    stop("countimp: the base method '", basis, "' for this 2lonly method was ",
         "not found.", call. = FALSE)
  impc <- f(y = yc, ry = ryc, x = xc, wy = !ryc)

  ## hand back to the rows
  voll <- yc
  voll[!ryc] <- as.numeric(impc)
  out <- voll[match(g, cl)]
  as.numeric(out[wy])
}


## Refuse row-wise imputation of a cluster-level variable ---------------------
##
## Called once per countimp() run, after the setup and before the first draw.
## A variable that is constant within a cluster and imputed row by row produces
## a completed data set the design cannot contain -- and, as measured at the top
## of this file, it does so without any visible sign.
##
## This is a behaviour change: a call that ran before now stops. That is
## deliberate and follows B73 (an observed zero under a zero-truncated method is
## an error, not a warning): the run would otherwise produce numbers that look
## right. The escape hatch is a named option rather than an argument, so it
## cannot be set by accident.
##
## Two guards against false alarms:
##   * only variables that are ACTUALLY imputed are checked (method non-empty),
##   * a variable with fewer than three observed clusters is skipped: with one
##     or two clusters, constancy is not evidence of anything.
.countimp_check_levels <- function(data, method, predictorMatrix,
                                   min_cluster = 3L) {
  if (isFALSE(getOption("countimp.check.levels", TRUE))) return(invisible(NULL))
  if (is.null(predictorMatrix) || !length(predictorMatrix)) return(invisible(NULL))
  vn <- colnames(predictorMatrix)
  ## `method` may be unnamed -- that is the common form in mice, and the one
  ## its own documentation uses: method = c("pmm", "", "", ""). The loop below
  ## indexes it BY NAME (method[[v]]), which then fails with "subscript out of
  ## bounds" before the check ever runs. An unnamed vector is positional with
  ## respect to vn; a named one is reordered onto vn and filled up, so that a
  ## partial vector -- method = c(y = "2l.nb2"), which mice allows -- selects
  ## the right positions. A length mismatch is a caller error and is reported
  ## as one, rather than silently checking the wrong variables.
  if (is.null(names(method))) {
    if (length(method) != length(vn))
      stop(sprintf(paste0(
        "countimp: length(method) = %d, but the predictorMatrix has %d ",
        "columns.\n  An unnamed method vector is read positionally, so the ",
        "two must have the same\n  length -- or name it: method = c(%s = ",
        "\"...\")."),
        length(method), length(vn), vn[1L]), call. = FALSE)
    names(method) <- vn
  } else {
    fehlt <- setdiff(vn, names(method))
    method <- c(method, stats::setNames(rep("", length(fehlt)), fehlt))[vn]
  }
  grp <- vn[apply(predictorMatrix == -2L, 2L, any)]
  if (!length(grp)) return(invisible(NULL))

  ## 2lonly methods are the correct treatment, so they are what we point AT,
  ## not what we complain about.
  for (v in vn[nzchar(method)]) {
    if (v %in% grp || grepl("^2lonly", method[[v]])) next
    y <- data[[v]]
    if (is.null(y) || !anyNA(y)) next
    for (g in grp) {
      gv <- data[[g]]
      if (is.null(gv)) next
      beob <- !is.na(y)
      if (length(unique(gv[beob])) < min_cluster) next
      if (!.countimp_constant_within(y, gv)) next
      stop("countimp: '", v, "' is constant within '", g,
           "', so it is a variable of that level,\n  but method \"",
           method[[v]], "\" imputes it row by row. The completed data would ",
           "then contain\n  clusters whose '", v, "' differs between their own ",
           "rows.\n",
           "  Use a cluster-level method instead, e.g. method = c(", v,
           " = \"2lonly.pmm\")\n  (also 2lonly.norm, 2lonly.pois, 2lonly.nb), ",
           "and code '", g, "' as -2 in that row\n  of the predictorMatrix. ",
           "To impute row-wise anyway: ",
           "options(countimp.check.levels = FALSE).",
           call. = FALSE)
    }
  }

  ## The mirror image: a 2lonly method pointed at too FINE a level.
  ##
  ## A school budget is constant within its school, and therefore within every
  ## class of that school as well. Aggregating to the class level passes every
  ## check the aggregation itself can make -- one value per class, all of them
  ## consistent -- and then draws one value PER CLASS. The completed data has a
  ## school with several budgets, and nothing in the output says so. Measured
  ## on 12 schools with 4 classes each: 4 of 12 schools came out with more than
  ## one budget, with no error and no warning.
  ##
  ## The core has refused this since 25 August (coarser_level_available() in
  ## clusterlevel.hpp); this is the R side catching up. It belongs HERE and not
  ## in .countimp_2lonly(), because a mice.impute.* function only ever sees the
  ## columns the predictorMatrix lets through -- measured: with the school
  ## coded 0 it is not in `x` at all, so the check could not see the level it
  ## has to warn about. This function gets the whole data frame.
  for (v in vn[nzchar(method)]) {
    if (!grepl("^2lonly", method[[v]])) next
    y <- data[[v]]
    if (is.null(y) || !anyNA(y)) next
    beob <- !is.na(y)
    ## the level the method will pick: the outermost -2 column the variable is
    ## constant within, which is what .countimp_2lonly_grp() does
    passend <- Filter(function(g) {
      gv <- data[[g]]
      !is.null(gv) && .countimp_constant_within(y, gv, beob)
    }, grp)
    if (!length(passend)) next
    gewaehlt <- passend[[which.min(vapply(passend,
      function(g) length(unique(data[[g]])), integer(1)))]]
    gw <- data[[gewaehlt]]

    for (k in setdiff(names(data), c(v, gewaehlt))) {
      kv <- data[[k]]
      if (is.null(kv) || length(kv) != length(gw)) next
      ## strictly coarser: the chosen grouping nested inside it, and fewer
      ## distinct values -- the same two conditions the core tests, and they
      ## are NOT redundant (a mere renaming is nested but not coarser)
      if (length(unique(kv)) >= length(unique(gw))) next
      if (!.countimp_is_nested(gw, kv)) next
      if (!.countimp_constant_within(y, kv, beob)) next
      stop("countimp: '", v, "' is constant on a COARSER level than the ",
           "grouping given: '", k, "'
  groups it as well, and '", gewaehlt,
           "' is nested inside it. Aggregating to the finer
  level draws one ",
           "value per inner group, so a single '", k, "' would end up with
  ",
           "several different values of '", v, "' -- a school with several ",
           "budgets.
  Code '", k, "' as -2 in that row of the ",
           "predictorMatrix instead.
",
           "  To group by the finer level anyway: ",
           "options(countimp.check.levels = FALSE).",
           call. = FALSE)
    }
  }
  invisible(NULL)
}
