## B31 -- the two-level hurdle family after consolidation.
##
## Up to 3.0.0 the 16 exported names in mice.impute.2l.hnb.R were 16 copies of
## one ~118-line body. Two defects were found in exactly that structure:
##   B30  the ZI family fitted a hurdle and drew untruncated counts
##   B31  the 8 .boot copies called ci(".countimp_rate")() without grp/obs_levels,
##        so units in clusters left out by the cluster bootstrap never received
##        u* ~ N(0, tau^2) -- glmmTMB returns a FINITE population-average value
##        for an unseen level, so the !is.finite(p) fallback never triggered.
##
## Both are copy-divergence defects: the fix was applied to some copies and not
## to others. The family now runs through ci(".countimp_2l_hurdle")(), so these tests
## guard the properties that the copies used to violate.
##
## Reference measurement: skripte/k01_referenz_hurdle.R. All 8 Bayesian
## variants are bit-identical across the rebuild (z = 0.000); the 4 hnb.boot
## variants return to the reference after the theta fix (max |z| = 0.95).

testthat::test_that("all 16 hurdle names exist and route through one engine", {
  nms <- c(outer(c("mice.impute.2l.hnb", "mice.impute.2l.hp"),
                 c("", ".noint.count", ".noint.zero", ".noint.both"), paste0))
  nms <- c(outer(nms, c("", ".boot"), paste0))
  testthat::expect_length(nms, 16L)
  for (n in nms) {
    testthat::expect_true(exists(n, mode = "function"), info = n)
    bod <- paste(deparse(body(get(n))), collapse = " ")
    ## every wrapper must delegate; none may carry its own model code
    testthat::expect_match(bod, ".countimp_2l_hurdle", fixed = TRUE, info = n)
    testthat::expect_false(grepl("glmmTMB::glmmTMB", bod, fixed = TRUE), info = n)
  }
})

testthat::test_that("the family and draw arguments reach the engine", {
  ## hp -> poisson, hnb -> nbinom2; .boot -> boot, else bayes
  chk <- function(n, fam, drw) {
    bod <- paste(deparse(body(get(n))), collapse = " ")
    testthat::expect_match(bod, paste0('"', fam, '"'), fixed = TRUE, info = n)
    testthat::expect_match(bod, paste0('"', drw, '"'), fixed = TRUE, info = n)
  }
  chk("mice.impute.2l.hnb",           "nbinom2", "bayes")
  chk("mice.impute.2l.hnb.boot",      "nbinom2", "boot")
  chk("mice.impute.2l.hp",            "poisson", "bayes")
  chk("mice.impute.2l.hp.noint.both.boot", "poisson", "boot")
})

testthat::test_that("the intercept flags are set per variant, not shared (B27)", {
  ## The .noint variants differ ONLY in the defaults of intercept.c/.z. This is
  ## the defect corrected in 2.6.0: intercept.z was declared but never used.
  fm <- function(n) formals(get(n))
  testthat::expect_true(fm("mice.impute.2l.hnb")$intercept.c)
  testthat::expect_true(fm("mice.impute.2l.hnb")$intercept.z)
  testthat::expect_false(fm("mice.impute.2l.hnb.noint.count")$intercept.c)
  testthat::expect_true( fm("mice.impute.2l.hnb.noint.count")$intercept.z)
  testthat::expect_true( fm("mice.impute.2l.hnb.noint.zero")$intercept.c)
  testthat::expect_false(fm("mice.impute.2l.hnb.noint.zero")$intercept.z)
  testthat::expect_false(fm("mice.impute.2l.hnb.noint.both")$intercept.c)
  testthat::expect_false(fm("mice.impute.2l.hnb.noint.both")$intercept.z)
})

testthat::test_that("B31: both model parts get grp and obs_levels", {
  ## The regression guard. In the engine source, every prediction helper call
  ## must carry the cluster arguments -- that is what B31 was missing.
  src <- paste(deparse(ci(".countimp_2l_hurdle")), collapse = "\n")
  ## Both helpers count: the Bayesian branch predicts through
  ## .countimp_rate(), the bootstrap branch through .countimp_boot_rate(),
  ## which additionally recovers the u_j of clusters the resample left out.
  ## Whichever is called, the cluster arguments must be there and NAMED --
  ## this guard caught the boot_rate call sites passing them positionally.
  muster <- c(".countimp_rate(", ".countimp_boot_rate(")
  n.rate <- sum(vapply(muster, function(m)
    lengths(regmatches(src, gregexpr(m, src, fixed = TRUE))), integer(1)))
  testthat::expect_gt(n.rate, 0L)
  ## no bare call: every occurrence is followed by grp = within the next 200
  ## characters (the call spans up to three lines)
  for (m in muster) {
    pos <- gregexpr(m, src, fixed = TRUE)[[1]]
    for (p in as.integer(pos)) {
      if (p < 0L) next
      seg <- substr(src, p, p + 200L)
      testthat::expect_match(seg, "grp = ", fixed = TRUE)
      testthat::expect_match(seg, "obs_levels", fixed = TRUE)
    }
  }
})

testthat::test_that("B31: the bootstrap draws theta instead of taking sigma()", {
  ## A cluster resample gives one point estimate of theta per fit, so the
  ## bootstrap must add theta's own draw -- and ci(".countimp_draw_theta")() also
  ## enforces the lower bound that guards against the B09a explosions.
  src <- paste(deparse(ci(".countimp_2l_hurdle")), collapse = "\n")
  testthat::expect_match(src, ".countimp_draw_theta", fixed = TRUE)
  testthat::expect_false(grepl("glmmTMB::sigma", src, fixed = TRUE))
})

testthat::test_that("a unit past the gate is never imputed as zero", {
  ## This is what separates a hurdle from zero inflation, and the direction of
  ## defect B30. Checked on the engine's own draw, with the gate forced open.
  testthat::skip_on_cran()
  set.seed(4242L)
  J <- 25L; nj <- 20L; N <- J * nj
  g <- rep(seq_len(J), each = nj)
  x1 <- stats::rnorm(N)
  u  <- stats::rnorm(J, 0, sqrt(0.3))[g]
  mu <- exp(1.0 + 0.6 * x1 + u)
  ## true hurdle: gate open with prob 0.75, positive part zero-truncated
  pos <- stats::rbinom(N, 1L, 0.75) == 1L
  y <- integer(N)
  z <- stats::rnbinom(sum(pos), size = 2, mu = mu[pos])
  while (any(z == 0L)) {
    j <- which(z == 0L); z[j] <- stats::rnbinom(length(j), size = 2, mu = mu[pos][j])
  }
  y[pos] <- z
  ry <- stats::rbinom(N, 1L, 0.75) == 1L
  x <- data.frame(g = g, x1 = x1)
  type <- c(-2L, 2L); names(type) <- colnames(x)

  imp <- suppressWarnings(
    mice.impute.2l.hnb(y = y, ry = ry, x = x, type = type))
  testthat::expect_length(imp, sum(!ry))
  testthat::expect_true(all(is.finite(imp)))
  testthat::expect_true(all(imp >= 0))
  testthat::expect_equal(imp, round(imp))          # counts, not continuous
  testthat::expect_gt(stats::sd(imp), 0)           # draws vary
  ## the imputations must contain BOTH zeros (gate closed) and positive counts
  testthat::expect_gt(mean(imp == 0), 0)
  testthat::expect_gt(mean(imp > 0), 0)
})

testthat::test_that("the engine honours type codes 3-6 per model part (B27)", {
  ## Codes 3/4 restrict a predictor to the count part, 5/6 to the zero part.
  ## The 2.6.0 defect was that the count part silently used ALL predictors.
  dec <- ci(".countimp_decode_type")(c(-2L, 3L, 5L, 1L), c("g", "only.c", "only.z", "both"))
  testthat::expect_equal(dec$group, "g")
  testthat::expect_true("only.c" %in% dec$c.fixed)
  testthat::expect_false("only.c" %in% dec$z.fixed)
  testthat::expect_true("only.z" %in% dec$z.fixed)
  testthat::expect_false("only.z" %in% dec$c.fixed)
  testthat::expect_true(all(c("both") %in% dec$c.fixed))
  testthat::expect_true(all(c("both") %in% dec$z.fixed))
})
