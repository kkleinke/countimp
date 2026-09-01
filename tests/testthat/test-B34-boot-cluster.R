## B34 -- a cluster bootstrap must draw the effect of clusters it did not
## sample.
##
## Sampling clusters with replacement means the bootstrap fit typically does
## not contain every cluster present in the target rows. For those glmmTMB with
## allow.new.levels = TRUE returns the population value, i.e. u_j = 0, instead
## of a draw from N(0, tau^2). Imputations for those units then vary too little
## between clusters -- which biases any downstream estimate of the random
## effect variance downwards, exactly the quantity multilevel data is about.
##
## Measured before/after the 3.0.0 consolidation (analyse/k02_referenz_2lcount_
## VORHER.csv vs _NACHHER.csv, 60 replications per method): the eight Bayesian
## variants are bit-identical (max |z| = 0.0000), all four bootstrap variants
## move up together and their distance to the Bayesian counterpart falls --
## 2l.nb.boot from 0.236 to 0.003. Bayes and bootstrap estimate the same thing.

testthat::test_that("B34: no draw path calls predict() directly", {
  ## The engine must route every prediction through ci(".countimp_predict_2l")(),
  ## which is where the unsampled-cluster draw lives.
  cand <- c("../../R/impute2lcount.R", "R/impute2lcount.R",
            "paket/R/impute2lcount.R", "../../../R/impute2lcount.R")
  f <- cand[file.exists(cand)][1]
  testthat::skip_if(is.na(f), "engine source not locatable")
  src <- sub("#.*$", "", readLines(f, warn = FALSE))
  src <- gsub('"[^"]*"', '""', src)
  direct <- grep("(^|[^._[:alnum:]])predict[[:space:]]*\\(", src, value = TRUE)
  ## ci(".countimp_predict_2l")() contains "predict" but is not a direct call
  direct <- direct[!grepl("countimp_predict", direct)]
  testthat::expect_equal(direct, character(0),
                         info = paste(direct, collapse = "\n"))
})

testthat::test_that("B34: every prediction carries the cluster arguments", {
  cand <- c("../../R/impute2lcount.R", "R/impute2lcount.R",
            "paket/R/impute2lcount.R", "../../../R/impute2lcount.R")
  f <- cand[file.exists(cand)][1]
  testthat::skip_if(is.na(f), "engine source not locatable")
  src <- paste(sub("#.*$", "", readLines(f, warn = FALSE)), collapse = " ")
  ## The engine predicts through two helpers, both of which take the cluster
  ## arguments: ci(".countimp_draw_2l")() on the Bayesian path and ci(".countimp_rate")()
  ## on the bootstrap path. Both ultimately call ci(".countimp_predict_2l")().
  helpers <- "\\.countimp_(draw_2l|rate|predict_2l)[[:space:]]*\\("
  sites <- gregexpr(helpers, src)[[1]]
  testthat::expect_true(sites[1] > 0)            # at least one call site
  n.sites <- length(sites)
  ## each call site must have both arguments before its closing paren; the
  ## engine is small enough that checking the whole text is sufficient
  n.grp <- length(gregexpr("grp[[:space:]]*=", src)[[1]])
  n.obs <- length(gregexpr("obs_levels[[:space:]]*=", src)[[1]])
  testthat::expect_gte(n.grp, n.sites)
  testthat::expect_gte(n.obs, n.sites)
})

testthat::test_that("B34: unsampled clusters get a drawn effect, not zero", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("glmmTMB")
  ## Fit on a subset of clusters, predict for a cluster the fit never saw.
  ## Repeated calls must differ -- that is the draw. If the code fell back to
  ## the population value they would all be identical.
  set.seed(4321)
  J <- 24L; nj <- 12L; N <- J * nj
  g <- rep(seq_len(J), each = nj)
  u <- rnorm(J, 0, 0.8)[g]
  x1 <- rnorm(N)
  y <- rpois(N, exp(1 + 0.4 * x1 + u))
  d <- data.frame(y = y, x1 = x1, g = factor(g))
  keep <- d$g %in% factor(1:20)                    # clusters 21-24 unseen
  fit <- suppressWarnings(glmmTMB::glmmTMB(y ~ 1 + x1 + (1 | g),
                                           data = d[keep, ], family = "poisson"))
  nd <- d[!keep, ][1:6, ]
  obs.levels <- unique(as.character(d$g[keep]))
  p <- replicate(12, as.numeric(ci(".countimp_predict_2l")(
         fit, nd, type = "link", grp = "g", obs_levels = obs.levels))[1])
  testthat::expect_gt(stats::sd(p), 0)             # it is a draw
  ## and the spread must be of the order of tau, not numerical noise
  testthat::expect_gt(stats::sd(p), 0.05)
})

testthat::test_that("B34: bootstrap and Bayes agree on the mean", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("glmmTMB")
  ## The two draw paths estimate the same quantity. Before the fix the
  ## bootstrap means sat systematically below their Bayesian counterparts;
  ## afterwards they agree to within Monte Carlo error. This is the property,
  ## not the numbers: a small run, generously bounded.
  set.seed(99)
  J <- 20L; nj <- 14L; N <- J * nj
  g <- rep(seq_len(J), each = nj)
  u <- rnorm(J, 0, 0.6)[g]
  x1 <- rnorm(N)
  y <- rpois(N, exp(1 + 0.4 * x1 + u))
  mis <- rbinom(N, 1L, 0.25) == 1L
  x <- data.frame(g = g, x1 = x1)
  ty <- c(-2L, 2L); names(ty) <- colnames(x)
  m.bay <- m.boo <- numeric(6)
  for (i in 1:6) {
    m.bay[i] <- mean(suppressWarnings(
      mice.impute.2l.poisson(y = y, ry = !mis, x = x, type = ty)))
    m.boo[i] <- mean(suppressWarnings(
      mice.impute.2l.poisson.boot(y = y, ry = !mis, x = x, type = ty)))
  }
  se <- sqrt(stats::var(m.bay) / 6 + stats::var(m.boo) / 6)
  testthat::expect_lt(abs(mean(m.bay) - mean(m.boo)), max(4 * se, 0.9))
})


## A left-out cluster is not an unseen cluster -------------------------------
##
## The block above fixed one half: a cluster the resample missed got u_j = 0
## instead of a draw. Reported from the simulation side on 31 August 2026, and
## reproduced here: the draw it got was from the PRIOR N(0, tau), which for
## that cluster is far too wide -- measured 0.6426 against a median posterior
## SD of 0.1203, a factor of 5.3. And its rows are not few: in the design of
## study 3, 18.1 of 50 clusters are missing from the average resample, which is
## 36.2 % of the rows to be imputed.
##
## The cluster is not unseen. Its observed rows are in the data; they merely
## missed the draw. The C++ core states the distinction and carries it in a
## parameter of its own -- `cluster_has_data` (glmm.hpp:496, branch at :526):
## a cluster WITH data gets p(u_j | y_j), only a cluster with no observed value
## at all gets N(0, sigma). The core is the reference, and these blocks hold
## the R side to it.
##
## Measured before and after, 40 clusters, m = 10, against the true u_j
## (scratch design of 31 August):
##
##             correlation with true u0   between-cluster SD
##   full data           0.9593                 0.3970
##   2l.nb2 (Bayes)      0.9467                 0.4054
##   2l.nb2.boot BEFORE  0.8803                 0.2444   <- 38 % too small
##   2l.nb2.boot AFTER   0.9492                 0.4202

testthat::test_that("B34: a left-out cluster keeps its own posterior", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("glmmTMB")
  ## The property, stated so it cannot pass by accident: imputations must
  ## carry the between-cluster spread of the data. The old code diluted it by
  ## replacing a third of the clusters' effects with prior draws.
  set.seed(51)
  J <- 30L; nj <- 40L; N <- J * nj
  g <- rep(seq_len(J), each = nj)
  u <- rnorm(J, 0, 0.7)
  x1 <- rnorm(N)
  y <- rpois(N, exp(1 + 0.4 * x1 + u[g]))
  mis <- rbinom(N, 1L, 0.3) == 1L
  x <- data.frame(g = g, x1 = x1)
  ty <- c(-2L, 2L); names(ty) <- colnames(x)

  ## cluster means of the imputations, averaged over draws, against the truth
  draw_ <- function(f) {
    M <- replicate(6, suppressWarnings(f(y = y, ry = !mis, x = x, type = ty)))
    rowMeans(apply(M, 2L, function(v) tapply(log1p(v), g[mis], mean)))
  }
  set.seed(7); kb <- draw_(mice.impute.2l.poisson)
  set.seed(7); kt <- draw_(mice.impute.2l.poisson.boot)
  wahr <- u[as.integer(names(kb))]
  ## The thresholds are measured, not guessed. Over six seeds in exactly this
  ## design, before and after the correction:
  ##
  ##                    before          after         threshold
  ##   correlation   0.799-0.873    0.916-0.976         0.90
  ##   sd ratio      0.674-0.777    0.935-1.030         0.86
  ##
  ## Both sit in the gap, with room on either side for Monte Carlo variation.
  ## An earlier version of this block used 0.85 and 0.6 -- picked by eye, and
  ## the mutation probe showed them to be blind: the pre-correction values
  ## passed them.
  testthat::expect_gt(stats::cor(kt, wahr), 0.90)
  ## and carry the between-cluster spread of the Bayesian branch, which draws
  ## u_j for every cluster and is the reference here
  testthat::expect_gt(stats::sd(kt), 0.86 * stats::sd(kb))
})

testthat::test_that("B34: a cluster with no observed value still draws from the prior", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("glmmTMB")
  ## The guard against over-correcting. Reading obs.levels from the full data
  ## must not turn a genuinely unobserved cluster into a known one: it has no
  ## conditional posterior, so N(0, tau) is right for it -- point 3 of the
  ## core, "a cluster nobody saw is not an average cluster, it is an unknown
  ## one". Repeated calls for such a cluster must therefore still differ.
  set.seed(88)
  J <- 22L; nj <- 20L; N <- J * nj
  g <- rep(seq_len(J), each = nj)
  u <- rnorm(J, 0, 0.8)
  x1 <- rnorm(N)
  y <- rpois(N, exp(1 + 0.4 * x1 + u[g]))
  mis <- g > 20L                       # clusters 21 and 22: nothing observed
  x <- data.frame(g = g, x1 = x1)
  ty <- c(-2L, 2L); names(ty) <- colnames(x)
  m <- replicate(10, mean(suppressWarnings(
        mice.impute.2l.poisson.boot(y = y, ry = !mis, x = x, type = ty))))
  testthat::expect_gt(stats::sd(m), 0)
})

testthat::test_that("B34: obs.levels is read from the full data, not the resample", {
  ## Source guard, in all four engine files. Reading it after the resample is
  ## precisely the defect: it makes .countimp_predict_2l() treat a left-out
  ## cluster as new. The order matters, so the test reads the order.
  dateien <- c("impute2lcount.R", "impute2lhurdle.R", "impute2lzi.R",
               "impute2lcmp.R")
  wurzel <- c("../../R", "R", "paket/R", "../../../R")
  wurzel <- wurzel[dir.exists(wurzel)][1]
  testthat::skip_if(is.na(wurzel), "engine source not locatable")
  for (dn in dateien) {
    src <- readLines(file.path(wurzel, dn), warn = FALSE)
    i.boot <- grep("<- .countimp_boot_clusters(", src, fixed = TRUE)
    i.lev  <- grep("obs.levels <- .countimp_obs_levels(", src, fixed = TRUE)
    testthat::expect_length(i.boot, 1L)
    testthat::expect_length(i.lev, 1L)
    ## the levels must come from dat.voll, and dat.voll must be taken before
    ## the resample overwrites dat
    testthat::expect_match(src[i.lev], "dat.voll", fixed = TRUE, info = dn)
    i.voll <- grep("^\\s*dat\\.voll <- dat\\s*$", src)
    testthat::expect_length(i.voll, 1L)
    testthat::expect_lt(i.voll, i.boot)
  }
})

testthat::test_that("B34: the refit holds every fixed parameter of the resample fit", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("glmmTMB")
  ## What makes the correction cheap and correct: the refit solves for u only.
  ## If any fixed parameter moved, the u_j would belong to a different model
  ## than the beta the bootstrap drew.
  set.seed(12)
  J <- 25L; nj <- 30L; N <- J * nj
  g <- rep(seq_len(J), each = nj)
  d <- data.frame(Y = rpois(N, exp(1 + rnorm(J, 0, .6)[g])),
                  x1 = rnorm(N), g = factor(g))
  form <- Y ~ 1 + x1 + (1 | g)
  rs <- ci(".countimp_boot_clusters")(d, "g")
  fit <- suppressWarnings(glmmTMB::glmmTMB(form, data = rs, family = "poisson"))
  rf <- ci(".countimp_boot_refit")(fit, form, d, ci(".countimp_family")("poisson", "log"))
  testthat::expect_true(rf$ok)
  testthat::expect_equal(unname(glmmTMB::fixef(rf$fit)$cond),
                         unname(glmmTMB::fixef(fit)$cond), tolerance = 1e-8)
  ## the refit knows every cluster, the resample fit does not
  testthat::expect_gt(length(rf$z$b), length(fit$obj$env$parList()$b))
  testthat::expect_equal(length(rf$z$b), J)
  ## and z is a draw, not a constant
  rf2 <- ci(".countimp_boot_refit")(fit, form, d, ci(".countimp_family")("poisson", "log"))
  testthat::expect_false(isTRUE(all.equal(rf$z$b, rf2$z$b)))
})

testthat::test_that("B34: the u_j are drawn, not taken as BLUPs", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("glmmTMB")
  ## The second half of the correction, and the one a between-cluster
  ## statistic cannot see: for clusters the resample DID sample, the old code
  ## took the BLUP -- a point estimate. The posterior spread around it
  ## (measured 0.1603 in the design of study 3) then appears in none of the m
  ## imputations, and that is the spread Rubin's rules need.
  ##
  ## Tested directly rather than through the imputations, because it is a
  ## sharp property: with the fit, the data and newdata held fixed, every step
  ## of .countimp_boot_rate() is deterministic EXCEPT the draw of u. So if
  ## two calls under different seeds agree, the draw is not happening. (The
  ## prior route for unseen clusters is excluded by keeping newdata inside the
  ## clusters the data has.)
  set.seed(19)
  J <- 20L; nj <- 30L; N <- J * nj
  g <- rep(seq_len(J), each = nj)
  d <- data.frame(Y = rpois(N, exp(1 + rnorm(J, 0, .6)[g])),
                  x1 = rnorm(N), g = factor(g))
  form <- Y ~ 1 + x1 + (1 | g)
  rs <- ci(".countimp_boot_clusters")(d, "g")
  fit <- suppressWarnings(glmmTMB::glmmTMB(form, data = rs, family = "poisson"))
  nd <- d[seq(1L, N, by = 7L), ]                  # every cluster is in the data
  fam <- ci(".countimp_family")("poisson", "log")
  ol <- ci(".countimp_obs_levels")(d, "g")
  r <- function(seed) {
    set.seed(seed)
    ci(".countimp_boot_rate")(fit, form, d, fam, nd, grp = "g", obs_levels = ol)
  }
  a <- r(1L); b <- r(2L)
  testthat::expect_length(a, nrow(nd))
  testthat::expect_true(all(is.finite(a)))
  testthat::expect_false(isTRUE(all.equal(a, b)))
  ## and the difference must be of the order of the posterior spread, not
  ## floating-point noise
  testthat::expect_gt(stats::sd(log(a) - log(b)), 0.01)
})

testthat::test_that("B34: the bootstrap branch falls back when the refit fails", {
  testthat::skip_on_cran()
  ## The correction must not turn a working imputation into an error. A
  ## non-glmmTMB fit is the cheapest way to make the refit refuse.
  rf <- ci(".countimp_boot_refit")(structure(list(), class = "lm"),
                                   y ~ 1, data.frame(y = 1), "poisson")
  testthat::expect_false(rf$ok)
  testthat::expect_match(rf$why, "glmmTMB")
})


test_that("B34: every bootstrap method actually gets its refit", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("glmmTMB")
  ## The guard this file was missing, and it cost a real defect: the fallback
  ## in .countimp_boot_rate() is deliberately quiet -- a failed refit must not
  ## turn a working imputation into an error -- and quiet is exactly how a
  ## half-applied correction survives. Measured on 31 August 2026: the hurdle
  ## methods ran ONE of their two refits. dat.voll was taken before `nz`, the
  ## response of the gate formula, was computed, so the gate refit failed for
  ## want of a column and degraded into the fallback without a sign. The
  ## imputations changed (so a before/after comparison looked convincing) --
  ## they were simply only half corrected.
  ##
  ## So: not "the numbers moved", but "the refit came through", for every
  ## method that has one. The diagnostics log is the channel, because that is
  ## where the fallback records itself.
  countimp_diagnostics(enable = TRUE, reset = TRUE)
  on.exit(countimp_diagnostics(enable = FALSE, reset = TRUE), add = TRUE)

  set.seed(77)
  J <- 25L; nj <- 24L; N <- J * nj
  g <- rep(seq_len(J), each = nj)
  u <- rnorm(J, 0, 0.6); x1 <- rnorm(N); z1 <- rnorm(N)
  y <- MASS::rnegbin(N, mu = exp(0.8 + u[g] + 0.3 * x1), theta = 2)
  y[runif(N) < 0.25] <- 0L                     # zeros for the two-part models
  mis <- rbinom(N, 1L, 0.28) == 1L
  x <- data.frame(g = g, x1 = x1, z1 = z1)
  ty <- c(-2L, 2L, 1L); names(ty) <- colnames(x)

  boot <- c("2l.poisson.boot", "2l.nb2.boot", "2l.cmp.boot", "2l.hp.boot",
            "2l.hnb.boot", "2l.zip.boot", "2l.zinb.boot")
  ## Fetched into a variable first, not `ci(".countimp_state")$x <- FALSE`:
  ## that is an assignment into the result of a function call and R refuses it
  ## with "target of assignment expands to non-language object" -- the same
  ## slip that silently killed two blocks before 3.0.0. An environment has
  ## reference semantics, so writing through the variable does reach it.
  zustand <- ci(".countimp_state")
  for (m in boot) {
    ## the note is written once per session, so clear the latch per method
    zustand$boot_u_noted <- FALSE
    countimp_diagnostics(enable = TRUE, reset = TRUE)
    set.seed(2024)
    v <- suppressWarnings(get(paste0("mice.impute.", m))(
           y = y, ry = !mis, x = x, type = ty))
    testthat::expect_true(all(is.finite(v)), info = m)
    dg <- countimp_diagnostics()
    faul <- if (is.null(dg) || !nrow(dg)) character(0) else
      dg$problems[grepl("boot_u_unavailable", dg$problems, fixed = TRUE)]
    testthat::expect_length(faul, 0L)
    if (length(faul)) testthat::fail(paste(m, "fell back:", faul[1L]))
  }
})
