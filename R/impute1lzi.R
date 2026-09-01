## One-level zero-inflated and hurdle count imputation -- one engine.
##
## The eight exported methods (zip, zinb, hp, hnb, each with a .boot variant)
## differ in exactly three switches: the pscl fitting function (zeroinfl or
## hurdle), the count distribution (poisson or negbin), and how the regression
## coefficients are obtained (posterior draw or bootstrap refit). Everything
## else -- building the two-part formula from `type`, assembling the data,
## predicting, drawing the imputations -- was copied eight times. This file
## holds it once; R/mice.impute.zip.R keeps the eight exported wrappers.
##
## Deliberately NOT changed here: the count part conditions on the point
## estimate of theta rather than drawing it. summary(fit)$vcov excludes
## log(theta), so the draw would have to come from fit$SE.logtheta separately.
## Measured in analyse/k03_theta_zi.csv (n = 400: SD ratios drawn/fixed 0.946
## and 0.978, both CIs covering 1) and analyse/k03b_theta_klein.csv
## (n = 60 ... 800, 500 replications each: ratios 0.987-1.412, all CIs but the
## one at n = 100 covering 1). So the omission has no systematic effect on the
## between-imputation variance, while a naive lognormal draw produces absurd
## values when theta is weakly determined (re-measured at n = 60: 37 of 500
## draws above 100, largest 7.7e294). Adding the draw requires the same
## lower-tail guard .countimp_draw_theta() applies; recorded as B36.

## Two-part formula from mice's `type` vector.
##   1 = both parts, 2 = count part only, 3 = zero part only, 0 = unused.
## An empty side becomes an intercept-only side rather than an error.
.countimp_zi_formula <- function(nam, type) {
  both  <- which(type == 1L)
  cnt   <- which(type == 2L)
  zro   <- which(type == 3L)
  count <- sort(unique(c(both, cnt)))
  zero  <- sort(unique(c(both, zro)))
  lhs <- if (length(count) == 0L) "Y ~ 1" else
    paste("Y ~", paste(nam[count], collapse = " + "))
  rhs <- if (length(zero) == 0L) "| 1" else
    paste("|", paste(nam[zero], collapse = " + "))
  stats::as.formula(paste(lhs, rhs))
}

## Posterior draw of the regression coefficients, pscl parameterisation.
##
## coef(fit) and summary(fit)$vcov agree in length and order for both zeroinfl
## and hurdle (checked for poisson and negbin: 4 coefficients, 4x4 vcov --
## log(theta) is NOT in summary(fit)$vcov, unlike vcov(fit)). The count slot is
## filled first, the zero slot with the remainder; nc is taken BEFORE the first
## assignment, because the shipped code re-read length(fit$coefficients$count)
## after overwriting it -- harmless while the lengths happened to match, but it
## made the split depend on its own side effect.
.countimp_zi_draw_beta <- function(fit) {
  s <- summary(fit)
  V <- s$vcov
  beta <- stats::coef(fit)
  if (is.null(V) || length(beta) != ncol(V))
    stop(sprintf(paste0(
      "countimp: pscl reports %d coefficients but a %s covariance matrix. ",
      "countimp draws beta from chol(summary(fit)$vcov) and needs the two to ",
      "agree. Please report this together with your pscl version (%s)."),
      length(beta), paste(dim(V), collapse = "x"),
      tryCatch(as.character(utils::packageVersion("pscl")),
               error = function(e) "unknown")), call. = FALSE)
  rv <- try(t(chol(V)), silent = TRUE)
  if (inherits(rv, "try-error"))
    stop("countimp: the covariance matrix of the imputation model is not ",
         "positive definite, so beta cannot be drawn. This usually means the ",
         "model is over-parameterised for the observed cases -- drop a ",
         "predictor or collapse a sparse category.", call. = FALSE)
  b.star <- beta + rv %*% stats::rnorm(ncol(rv))
  nc <- length(fit$coefficients$count)
  nz <- length(fit$coefficients$zero)
  if (nc + nz != length(b.star))
    stop("countimp: the count and zero coefficient slots hold ", nc, " + ", nz,
         " values but the draw has ", length(b.star),
         ". This is a bug in countimp -- please report it.", call. = FALSE)
  ## keep the slot names: pscl's predict() indexes theta by name, and while it
  ## uses the coefficients positionally, preserving names keeps the fit
  ## printable and inspectable after the draw
  cc <- b.star[seq_len(nc)]
  zz <- b.star[nc + seq_len(nz)]
  names(cc) <- names(fit$coefficients$count)
  names(zz) <- names(fit$coefficients$zero)
  fit$coefficients$count <- cc
  fit$coefficients$zero  <- zz
  fit
}

#' One-level zero-inflated and hurdle imputation engine
#'
#' Internal engine behind \code{mice.impute.zip()}, \code{mice.impute.zinb()},
#' \code{mice.impute.hp()}, \code{mice.impute.hnb()} and their \code{.boot}
#' variants. Not intended to be called directly.
#'
#' @param y Incomplete numeric vector of counts.
#' @param ry Logical vector, \code{TRUE} where \code{y} is observed.
#' @param x Matrix or data frame of predictors, complete.
#' @param type Predictor codes: 1 both parts, 2 count part, 3 zero part.
#' @param wy Logical vector marking the cases to impute; defaults to
#'   \code{!ry}.
#' @param engine \code{"zeroinfl"} or \code{"hurdle"}.
#' @param dist \code{"poisson"} or \code{"negbin"} for the count part.
#' @param draw \code{"bayes"} for a posterior draw of the coefficients,
#'   \code{"boot"} for a bootstrap refit.
#' @return Numeric vector of imputations, one per \code{TRUE} in \code{wy}.
#' @keywords internal
.countimp_1l_zi <- function(y, ry, x, type, wy = NULL,
                            engine = c("zeroinfl", "hurdle"),
                            dist = c("poisson", "negbin"),
                            draw = c("bayes", "boot")) {
  engine <- match.arg(engine)
  dist   <- match.arg(dist)
  draw   <- match.arg(draw)
  if (is.null(wy)) wy <- !ry

  ## drop = FALSE throughout: with a single predictor x[ry, ] collapses to a
  ## vector, and the shipped code then relied on data.frame() naming the column
  ## "X" on both sides so the mismatch cancelled out. Keeping the frame keeps
  ## the real predictor names, which is what the formula and predict() need.
  X <- as.data.frame(x[ry, , drop = FALSE])
  nam <- colnames(X)
  dat <- data.frame(Y = y[ry], X)
  form <- .countimp_zi_formula(nam, type)

  ## predict() needs the predictor names the formula was built with
  newdata <- as.data.frame(x[wy, , drop = FALSE])
  colnames(newdata) <- nam

  ## Zero-free data checked HERE and not inside .countimp_zi_fit(): raised from
  ## inside the closure, .countimp_draw_retry() repeats it three times and then
  ## reports a convergence failure ("extreme predictor values, a separated
  ## model ...") with the actual cause buried in a parenthesis -- exactly the
  ## wrong diagnosis for a user who picked the wrong method. Same reasoning and
  ## same placement as .countimp_zt_check_zeros() for ztp/ztnb (B73).
  .countimp_zi_check_zeros(y[ry], engine)

  ## One draw, repeatable. Same contract as in .countimp_1l_count(): refit or
  ## redraw beta from the ambient RNG state on every call, and hand fit and mu
  ## to countimp_check() so a hard failure is retried rather than fatal (B56).
  ein_zug <- function() {
    d_i <- dat
    if (identical(draw, "boot")) {
      n.obs <- nrow(d_i)
      d_i <- d_i[sample.int(n.obs, n.obs, replace = TRUE), , drop = FALSE]
    }
    fit <- .countimp_zi_fit(form, d_i, engine, dist)
    if (identical(draw, "bayes")) fit <- .countimp_zi_draw_beta(fit)
    mu <- suppressWarnings(try(as.numeric(stats::predict(
      fit, newdata = newdata, type = "count", na.action = stats::na.pass)),
      silent = TRUE))
    if (inherits(mu, "try-error")) mu <- NULL
    list(imp = .countimp_rzi(fit, newdata), fit = fit, mu = mu)
  }

  .countimp_draw_retry(ein_zug, y_obs = y[ry],
                       method = paste(engine, dist, sep = "."))
}

## The single place where pscl's fitting functions are called. Everything the
## package knows about zeroinfl/hurdle argument names lives here, so a change
## on pscl's side is a one-line fix rather than eight.
## No zeros means no zero process to model, and pscl::zeroinfl()/hurdle() refuse
## outright with "invalid dependent variable, minimum count is not zero". The
## check is named and separate so the message can say WHICH method to use
## instead -- the mirror image of .countimp_zt_check_zeros() (B73).
.countimp_zi_check_zeros <- function(y, engine) {
  y <- y[!is.na(y)]
  if (!length(y) || any(y == 0)) return(invisible(TRUE))
  stop("countimp: the ", engine, " methods model a separate zero process, but ",
       "none of the ", length(y), " observed value(s) is zero.\n",
       "  Without zeros there is nothing for the zero part to fit. Use ",
       "\"ztp\"/\"ztnb\" (zero-truncated)\n",
       "  when zeros cannot occur, or \"poisson\"/\"nb\" when they are ordinary ",
       "counts.", call. = FALSE)
}


.countimp_zi_fit <- function(form, dat, engine, dist) {
  ## Check availability BEFORE the call. Without this, a missing pscl surfaces
  ## as "failed to fit the imputation model", which points the user at their
  ## data instead of at the missing package (B43).
  .countimp_need_pkg("pscl",
    what = "the single-level zero-inflated and hurdle methods",
    ersatz = "the multilevel variants (2l.zip, 2l.zinb, 2l.hp, 2l.hnb), which use glmmTMB")
  fit <- try(switch(engine,
    zeroinfl = pscl::zeroinfl(form, data = dat, dist = dist, link = "logit"),
    hurdle   = pscl::hurdle(form, data = dat, dist = dist, link = "logit")),
    silent = TRUE)
  if (inherits(fit, "try-error"))
    stop(sprintf(paste0(
      "countimp: pscl::%s() failed to fit the imputation model (%s count part). ",
      "The original message was: %s"), engine, dist,
      ## conditionMessage(), NOT sub() on as.character(): the pattern
      ## "^Error[^:]*:" stops at the FIRST colon, and a try-error's text carries
      ## the failing CALL before the message. Measured on pscl::zeroinfl()
      ## refusing zero-free data, the old form produced
      ## ":zeroinfl(y ~ x, dist = ..., data = ... :  invalid dependent ..."
      ## -- the call fragment survived and the reason was buried. The condition
      ## object carries the message alone: "invalid dependent variable,
      ## minimum count is not zero".
      trimws(conditionMessage(attr(fit, "condition")))), call. = FALSE)
  fit
}
