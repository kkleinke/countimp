## impute2lzi.R -- one engine for the two-level zero-inflated methods
##
## Up to 2.6.0 the 16 exported names in mice.impute.2l.zinb.R were 16 copies of
## the same 140-line body, differing in four places: the count family
## (nbinom2 / poisson), the draw (Bayes / cluster bootstrap), and the two
## intercept flags -- which are already function arguments, so the .noint
## variants only fix their defaults. That left 2338 lines in which a fix had to
## be applied 16 times, and two defects corrected in 2.6.0 were exactly
## fixes that had been applied to some copies and not others.
##
## THE DEFECT THIS FILE CORRECTS (B30). The old body fit a HURDLE and drew as
## if the count part were unrestricted:
##
##   stage 1  binomial GLMM on 1{Y = 0}          -> P(structural zero)
##   stage 2  nbinom2 GLMM on the POSITIVE subset -> mu
##   draw     rbinom(1 - p) to decide "positive", then rnbinom(mu) -- which
##            returns 0 with probability (theta/(theta+mu))^theta
##
## So a unit assigned to the count process could still be imputed as zero, and
## the count model was fitted to positive-only data without accounting for the
## truncation. Both parts are wrong in the same direction: too many zeros, too
## small a conditional mean.
##
## Measured (skripte/b30_zi_comparison.R, R = 200, m = 5, n = 625, 38-40 % MAR,
## coverage of the marginal x effect, nominal 95 %):
##
##   true process        current   hurdle throughout   genuine ZI
##   zero-inflated        0.885*         0.890*          0.915
##   hurdle               0.870*        0.920            0.870*
##   (* outside 2 Monte Carlo errors of 0.95)
##
## Each specification is best under its own process and the old code is worse
## under both. The methods are NAMED zip/zinb, and the hurdle methods exist
## separately as 2l.hp/2l.hnb, so these now fit what their name says: a genuine
## zero-inflation model, one glmmTMB fit with a ziformula, both linear
## predictors moved by ONE draw from the joint posterior
## (.countimp_draw_2l_zi()), and the count draw untruncated -- which is correct
## here, because under zero inflation a count-process unit may legitimately
## be zero.

## Is `innen` nested within `aussen`? True when every level of the inner factor
## occurs under exactly one level of the outer one -- classes inside schools.
## Crossed factors (pupils by teacher, say) fail this, and a cluster bootstrap
## is not defined for them: resampling one factor breaks the other's design.
## Is `inner` nested inside `outer`? Every inner unit must belong to exactly
## one outer unit -- classes inside schools, yes; classes inside teachers who
## teach across schools, no. Takes the two VECTORS, so that callers holding
## columns and callers holding names can both use it.
##
## The counterpart of is_nested_in() in the core's clusterlevel.hpp.
.countimp_is_nested <- function(inner, outer) {
  if (length(inner) != length(outer)) return(FALSE)
  tab <- table(as.character(inner), as.character(outer))
  all(rowSums(tab > 0) == 1L)
}


## Which grouping level is the outermost? Decided from the DATA, not from the
## order the formula happens to list the terms in: the outer factor is the one
## every other is nested within. Returns NA when the factors are not a nesting
## chain.
.countimp_outermost_level <- function(dat, grp) {
  if (length(grp) == 1L) return(grp)
  ist_aussen <- vapply(grp, function(a)
    all(vapply(setdiff(grp, a), function(i)
      .countimp_is_nested(dat[[i]], dat[[a]]), logical(1))), logical(1))
  if (sum(ist_aussen) != 1L) return(NA_character_)
  grp[which(ist_aussen)]
}


## Cluster bootstrap: resample whole clusters, not rows. Returns the resampled
## data frame. Kept separate so the Bayes and boot paths differ in one line.
##
## With nested levels the resample is HIERARCHICAL: draw whole schools, and
## every class of a drawn school comes with it. That is the only version that
## preserves the design -- resampling classes independently would break their
## nesting inside schools, and resampling schools while re-labelling their
## classes would invent clusters. Because the inner levels travel with their
## outer one, this is exactly "resample the outermost factor and take its rows",
## which is what the single-level code already did; the work is in FINDING the
## outermost factor and in refusing the case where there is none.
##
## A school drawn twice appears with the same identifier twice, as in the
## single-level version, so the two copies are fitted as one cluster. That is
## the established behaviour of this package's bootstrap (B34), not a decision
## taken here.
.countimp_boot_clusters <- function(dat, grp.name) {
  if (length(grp.name) > 1L) {
    aussen <- .countimp_outermost_level(dat, grp.name)
    if (is.na(aussen))
      stop("countimp: the cluster bootstrap needs the grouping factors to be ",
           "nested, but\n  ", paste(grp.name, collapse = " and "),
           " are crossed (no factor contains all the others).\n",
           "  Resampling one of them would break the other's design. Use the ",
           "Bayesian\n  variant (the method name without '.boot') for a ",
           "crossed design.", call. = FALSE)
    grp.name <- aussen
  }
  cls <- sample(unique(dat[[grp.name]]), replace = TRUE)
  idx <- unlist(lapply(cls, function(g) which(dat[[grp.name]] == g)),
                use.names = FALSE)
  out <- dat[idx, , drop = FALSE]
  rownames(out) <- NULL
  out
}


## The random effects belong to the DATA, not to the resample ----------------
##
## A cluster bootstrap leaves clusters out -- that is what sampling with
## replacement does, and it is not a defect. Measured in the design of study 3
## (50 clusters, sizes 100-400): 18.1 clusters are missing from the average
## resample, which is 36.2 % of the rows to be imputed.
##
## What was wrong is how those rows were then treated. `obs.levels` came from
## the RESAMPLE, so a left-out cluster counted as an unseen level and
## .countimp_predict_2l() drew its effect from the prior N(0, tau) -- measured
## 5.3 times wider than that cluster's own posterior (tau = 0.6426 against a
## median posterior SD of 0.1203). But the cluster is not unseen: its observed
## rows are right there in the data, they merely missed the draw.
##
## The C++ core states the distinction the R side had lost, and carries it in a
## parameter of its own (`cluster_has_data`, glmm.hpp:496, branch at :526):
## a cluster WITH data gets its conditional posterior p(u_j | y_j); only a
## cluster with no observed value at all gets N(0, sigma) -- "a cluster nobody
## saw is not an average cluster, it is an unknown one". The core is the
## reference here, so this function brings the R side to it.
##
## What the bootstrap keeps is its own job: (beta, Sigma, theta) come from the
## resample fit. The u_j are then solved on ALL observed data with those
## parameters held fixed, and drawn from their conditional posterior -- the
## same factorisation the Bayesian branch uses, only with a different source
## for the parameters. Two consequences worth naming:
##
##   * A cluster with no observed y is still absent from the fixed fit, so it
##     still gets the prior draw. That is point 3 of the core, and correct.
##   * Fixing every parameter makes the refit cheap: measured 0.25 s against
##     1.01 s for the full fit (factor 0.25), because it is the inner
##     optimisation step the full fit already solves in each iteration.
##
## Returns the rate on the response scale, i.e. it replaces .countimp_rate()
## in the bootstrap branch. On any failure it falls back to that function, so
## the branch degrades to its previous behaviour instead of aborting.
.countimp_boot_refit <- function(fit, form, data, family, ziformula = NULL) {
  if (!inherits(fit, "glmmTMB")) return(list(ok = FALSE, why = "fit is not a glmmTMB object"))
  pl <- tryCatch(fit$obj$env$parList(), error = function(e) NULL)
  if (is.null(pl) || !length(pl$b))
    return(list(ok = FALSE, why = "no random-effect block in the fit"))

  ## Freeze every FIXED parameter; the random-effect blocks stay free. A
  ## zero-inflated or dispersion-modelled fit has more than one of them (b,
  ## bzi, bdisp), and freezing those would defeat the purpose. `start` carries
  ## the resample estimates in, `map` keeps the optimiser from moving them.
  frei <- c("b", "bzi", "bdisp")
  fest <- setdiff(names(pl), frei)
  fest <- fest[lengths(pl[fest]) > 0L]
  karte <- lapply(pl[fest], function(v) factor(rep(NA, length(v))))
  args <- list(formula = form, data = data, family = family,
               start = pl[fest], map = karte)
  if (!is.null(ziformula)) args$ziformula <- ziformula
  ff <- tryCatch(suppressWarnings(do.call(glmmTMB::glmmTMB, args)),
                 error = function(e) NULL)
  if (is.null(ff)) return(list(ok = FALSE, why = "the fixed-parameter refit failed"))

  ## Draw the random effects from their conditional posterior. With every
  ## fixed parameter mapped out, the joint precision is the precision of the
  ## random-effect blocks alone.
  P <- suppressWarnings(tryCatch(
        TMB::sdreport(ff$obj, getJointPrecision = TRUE)$jointPrecision,
        error = function(e) NULL))
  if (is.null(P)) return(list(ok = FALSE, why = "sdreport(getJointPrecision) failed"))
  nmz <- colnames(as.matrix(P))
  benannt <- !is.null(nmz) && any(nzchar(nmz))
  pl2 <- ff$obj$env$parList()
  laengen <- lengths(pl2[frei])
  if (benannt) {
    idx <- lapply(frei, function(k) which(nmz == k))
  } else {
    ## unnamed: the blocks follow parList() order, b first
    if (nrow(P) != sum(laengen))
      return(list(ok = FALSE, why = "precision does not match the random-effect blocks"))
    ende <- cumsum(laengen); anf <- ende - laengen + 1L
    idx <- lapply(seq_along(frei), function(k)
      if (laengen[k]) seq.int(anf[k], ende[k]) else integer(0))
  }
  names(idx) <- frei
  if (length(idx$b) != laengen[["b"]])
    return(list(ok = FALSE, why = "precision does not match the b block"))

  alle <- unlist(idx, use.names = FALSE)
  U <- suppressWarnings(tryCatch(Matrix::chol(P[alle, alle, drop = FALSE]),
                                 error = function(e) NULL))
  if (is.null(U))
    return(list(ok = FALSE, why = "the posterior precision is not positive definite"))
  z <- tryCatch(as.numeric(Matrix::solve(U, stats::rnorm(length(alle)))),
                error = function(e) NULL)
  if (is.null(z) || !all(is.finite(z)))
    return(list(ok = FALSE, why = "solving the Cholesky system failed"))

  ## split z back into the blocks, in the order they were stacked
  pos <- cumsum(vapply(idx, length, integer(1)))
  teil <- list(); von <- 1L
  for (k in seq_along(idx)) {
    teil[[frei[k]]] <- if (length(idx[[k]])) z[von:pos[k]] else numeric(0)
    von <- pos[k] + 1L
  }
  list(ok = TRUE, fit = ff, z = teil)
}

## The rate for the bootstrap branch: replaces .countimp_rate(). On any
## failure it falls back to that function, so the branch degrades to its
## previous behaviour instead of aborting.
.countimp_boot_rate <- function(fit, form, data, family, newdata, grp,
                                obs_levels) {
  zurueck <- function(why) {
    if (!isTRUE(.countimp_state$boot_u_noted)) {
      .countimp_state$boot_u_noted <- TRUE
      .countimp_note_event("boot_u_unavailable", why)
    }
    .countimp_rate(fit, newdata, grp = grp, obs_levels = obs_levels)
  }
  rf <- .countimp_boot_refit(fit, form, data, family)
  if (!isTRUE(rf$ok)) return(zurueck(rf$why))

  ## eta from the refit -- which knows every observed cluster, so only truly
  ## unobserved ones take the prior route inside .countimp_predict_2l() -- and
  ## then the shift from BLUP to draw. Rows of an unobserved cluster have a
  ## zero row in Z, so the shift leaves their prior draw untouched.
  eta <- tryCatch(.countimp_predict_2l(rf$fit, newdata, type = "link",
                                       grp = grp, obs_levels = obs_levels),
                  error = function(e) NULL)
  if (is.null(eta)) return(zurueck("prediction from the refit failed"))
  Z <- .countimp_re_design(rf$fit, newdata)
  if (is.null(Z) || ncol(Z) != length(rf$z$b))
    return(zurueck("random-effect design for newdata could not be rebuilt"))
  eta <- as.numeric(eta) + as.numeric(Z %*% rf$z$b)

  linkinv <- tryCatch(stats::family(rf$fit)$linkinv, error = function(e) exp)
  if (!is.function(linkinv)) linkinv <- exp
  as.numeric(linkinv(eta))
}

## The engine. `family` is "nbinom2" (zinb) or "poisson" (zip); `draw` is
## "bayes" or "boot".
.countimp_2l_zi <- function(y, ry, x, type,
                            family = c("nbinom2", "poisson"),
                            draw = c("bayes", "boot"),
                            intercept.c = TRUE, intercept.z = TRUE,
                            wy = NULL, EV = FALSE) {
  family <- match.arg(family)
  draw   <- match.arg(draw)
  if (is.null(wy)) wy <- !ry

  X   <- data.frame(x[ry, , drop = FALSE])
  nam <- colnames(X)
  dec <- .countimp_decode_type(type, nam)
  dat <- data.frame(Y = y[ry], X)

  ## Clusters present in the DATA, so units from genuinely unseen clusters get
  ## u* ~ N(0, tau^2) instead of the population average (see predict2l.R).
  ## Read from the full data, not from the resample: a cluster the resample
  ## left out still has observed rows, and .countimp_boot_refit() recovers its
  ## conditional posterior from them.
  dat.voll <- dat
  if (identical(draw, "boot")) dat <- .countimp_boot_clusters(dat, dec$group)
  obs.levels <- .countimp_obs_levels(dat.voll, dec$group)

  newdata <- data.frame(x[wy, , drop = FALSE])
  colnames(newdata) <- nam

  ## One fit, both parts. The zero part carries its own random structure, so
  ## intercept.z and the type codes 5/6 act on it and nothing else -- the
  ## defect corrected in 2.6.0 was that intercept.c governed both.
  f.cond <- .countimp_2l_formula(dec, "count", response = "Y",
                                 intercept = intercept.c)
  f.zi   <- .countimp_2l_formula(dec, "zero",  response = "",
                                 intercept = intercept.z)
  fit <- glmmTMB::glmmTMB(formula = f.cond, ziformula = f.zi, data = dat,
                          family = .countimp_family(family, "log"))

  if (identical(draw, "bayes")) {
    dr <- .countimp_draw_2l_zi(fit, newdata, grp = dec$group,
                               obs_levels = obs.levels)
  } else {
    ## The bootstrap carries the uncertainty in the parameters; the random
    ## effects belong to the data and are drawn from their conditional
    ## posterior at those parameters -- for BOTH parts, which is why the refit
    ## keeps b and bzi free. Falls back to the plain prediction if the refit
    ## does not come through.
    rf <- .countimp_boot_refit(fit, f.cond, dat.voll,
                               .countimp_family(family, "log"), f.zi)
    if (!isTRUE(rf$ok)) {
      if (!isTRUE(.countimp_state$boot_u_noted)) {
        .countimp_state$boot_u_noted <- TRUE
        .countimp_note_event("boot_u_unavailable", rf$why)
      }
      fz <- fit; zc <- NULL; zz <- NULL
    } else {
      fz <- rf$fit; zc <- rf$z$b; zz <- rf$z$bzi
    }
    eta <- as.numeric(.countimp_predict_2l(fz, newdata, type = "link",
                                           grp = dec$group,
                                           obs_levels = obs.levels))
    Zc <- if (is.null(zc)) NULL else .countimp_re_design(fz, newdata, "cond")
    if (!is.null(Zc) && ncol(Zc) == length(zc)) eta <- eta + as.numeric(Zc %*% zc)
    etaz <- as.numeric(stats::predict(fz, newdata = newdata, type = "zlink",
                                      allow.new.levels = TRUE))
    Zz <- if (is.null(zz) || !length(zz)) NULL else
            .countimp_re_design(fz, newdata, "zi")
    if (!is.null(Zz) && ncol(Zz) == length(zz)) etaz <- etaz + as.numeric(Zz %*% zz)
    dr <- list(
      mu = exp(eta),
      pi = stats::plogis(etaz),
      theta = suppressWarnings(as.numeric(.countimp_theta(fit, pflicht = FALSE) %||% NA_real_)))
  }

  ## structural zero with probability pi; otherwise a count draw, which may
  ## itself be zero -- that is what distinguishes ZI from a hurdle
  structural <- stats::rbinom(length(dr$pi), size = 1L, prob = dr$pi) == 1L
  cnt <- if (identical(family, "poisson")) {
    stats::rpois(length(dr$mu), lambda = dr$mu)
  } else {
    stats::rnbinom(length(dr$mu), size = dr$theta, mu = dr$mu)
  }
  imp <- ifelse(structural, 0L, cnt)

  if (EV) imp <- .countimp_ev_handle(imp, y, x, wy)
  imp
}
