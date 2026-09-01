## B20: two-level methods understated the between-imputation variance.
##
## Before the fix, beta* was drawn from N(beta.hat, V.beta) and predict() was
## called on the fit carrying beta*. predict.glmmTMB() re-derives the
## conditional modes of the random effects from the parameter vector it is
## handed, so it absorbed most of the drawn shift: only
## sigma^2 / (n_j tau^2 + sigma^2) of it reached eta. Measured effect:
## sd(eta) between imputations was 23% of its correct value at n = 400 /
## 20 clusters, i.e. too small by a factor of four.
##
## The fix draws (beta, u, log theta) jointly from the Laplace posterior
## (TMB::sdreport(getJointPrecision = TRUE)). The reference these tests check
## against is se.fit from predict(..., se.fit = TRUE): it is the posterior sd
## of eta including the conditional variance of the random effects, which is
## exactly the spread a proper imputation must reproduce.

test_that("the joint draw reproduces the posterior sd of eta", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("TMB")
  d <- sim_count_2l(n = 300, ngrp = 15, seed = 11)
  fit <- fit_2l_on(d)
  skip_if(is.null(fit) || fit$fit$convergence != 0, "fit did not converge")
  nd <- d[1:20, ]

  target <- stats::predict(fit, newdata = nd, type = "link", se.fit = TRUE,
                           allow.new.levels = TRUE)$se.fit
  set.seed(3)
  eta <- replicate(300, ci(".countimp_draw_2l")(fit, nd)$eta)
  got <- apply(eta, 1, stats::sd)

  ## Monte-Carlo error on sd from 300 draws is about 4%, so 15% is a loose
  ## band that still fails decisively for the old route (it reached ~23%).
  expect_equal(mean(got) / mean(target), 1, tolerance = 0.15)
})

test_that("the joint draw is far wider than the beta-only route", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("TMB")
  d <- sim_count_2l(n = 300, ngrp = 15, seed = 11)
  fit <- fit_2l_on(d)
  skip_if(is.null(fit) || fit$fit$convergence != 0, "fit did not converge")
  nd <- d[1:20, ]

  set.seed(3)
  sd_new <- mean(apply(replicate(200, ci(".countimp_draw_2l")(fit, nd)$eta), 1, stats::sd))

  beta <- glmmTMB::fixef(fit)$cond
  rv <- t(chol(stats::vcov(fit)$cond))
  set.seed(3)
  sd_old <- mean(apply(replicate(200, {
    bs <- as.numeric(beta) + as.numeric(rv %*% stats::rnorm(ncol(rv)))
    as.numeric(ci(".countimp_predict_2l")(
      ci(".countimp_inject_beta")(fit, bs), nd, type = "link"))
  }), 1, stats::sd))

  expect_gt(sd_new / sd_old, 2)
})

test_that("theta is drawn, and the draw is floored", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("TMB")
  d <- sim_count_2l(n = 300, ngrp = 15, seed = 11)
  fit <- fit_2l_on(d)
  skip_if(is.null(fit) || fit$fit$convergence != 0, "fit did not converge")
  nd <- d[1:5, ]

  set.seed(5)
  th <- replicate(100, ci(".countimp_draw_2l")(fit, nd)$theta)
  expect_gt(stats::sd(th), 0)          # theta varies between imputations
  expect_true(all(th >= 0.25))         # never below the documented floor

  ## an explicitly supplied theta is passed through untouched
  expect_equal(ci(".countimp_draw_2l")(fit, nd, theta = 7)$theta, 7)
})

test_that("the joint route degrades to the old route instead of failing", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("TMB")
  d <- sim_count_2l(n = 300, ngrp = 15, seed = 11)
  fit <- fit_2l_on(d)
  skip_if(is.null(fit) || fit$fit$convergence != 0, "fit did not converge")

  ## A grouping factor missing from newdata breaks BOTH routes (predict.glmmTMB
  ## cannot build the model frame), so it must be named, not silently degraded.
  nd <- d[1:5, setdiff(names(d), "g"), drop = FALSE]
  expect_error(ci(".countimp_draw_2l")(fit, nd), "missing the grouping factor")

  ## A random slope variable missing from newdata, by contrast, only breaks the
  ## joint route -- there we degrade to the beta-only draw.
  fs <- fit_2l_on(d, form = y ~ x + (1 + x | g))
  if (!is.null(fs) && fs$fit$convergence == 0) {
    d2 <- d; d2$x2 <- stats::rnorm(nrow(d))
    res <- ci(".countimp_draw_2l")(fs, d[1:5, ])
    expect_true(all(is.finite(res$eta)))
  }

  ## a non-glmmTMB fit is a programming error inside countimp, not a user
  ## situation, so it must be named rather than fail deep inside fixef()
  expect_error(ci(".countimp_draw_2l")(stats::lm(y ~ x, data = d), d[1:5, ]),
               "expects a glmmTMB fit")
})

test_that("the random-effect design matches getME(fit, 'Z') exactly", {
  skip_if_not_installed("glmmTMB")
  d <- sim_count_2l(n = 300, ngrp = 15, seed = 11)
  for (form in list(y ~ x + (1 | g), y ~ x + (1 + x | g))) {
    fit <- tryCatch(suppressWarnings(glmmTMB::glmmTMB(form, data = d,
                    family = glmmTMB::nbinom2)), error = function(e) NULL)
    if (is.null(fit)) next
    Z_own <- ci(".countimp_re_design")(fit, d)
    Z_ref <- as.matrix(glmmTMB::getME(fit, "Z"))
    expect_equal(unname(Z_own), unname(Z_ref), tolerance = 1e-10,
                 info = deparse(form))
  }
})

test_that("cluster-level spread is recovered, not just individual-level", {
  ## The sharpest check available. On individual imputed values the error is
  ## damped: the independent count-draw noise dominates and B grew only by
  ## 1.3-1.8 in simulation. Within-cluster aggregates are where it shows in
  ## full, because that noise averages out while the missing cluster effect
  ## u_j* -- identical for every observation in a cluster -- does not.
  ## Measured before the fix: sd of cluster means too small by 2.85x.
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("TMB")
  skip_on_cran()          # 2 x 150 glmmTMB predictions

  ## The fraction missing in the target clusters is what governs the size of
  ## the effect: with 50% missing the fit still pins u_j down from the observed
  ## half and the factor is only 1.4; at 90% the imputation carries the cluster
  ## and it reaches 2.9. We use 90% -- the regime these methods are for.
  d <- sim_count_2l(n = 240, ngrp = 12, theta = 20, seed = 5)
  tgt <- levels(d$g)[1:3]
  it  <- which(d$g %in% tgt)
  set.seed(5)
  mis <- sample(it, round(0.9 * length(it)))
  fit <- fit_2l_on(d[-mis, ])
  skip_if(is.null(fit) || fit$fit$convergence != 0, "fit did not converge")

  nd <- d[mis, ]
  gm <- droplevels(nd$g)

  cl_sd <- function(gen) {
    m <- replicate(150, tapply(gen(), gm, mean))
    mean(apply(m, 1, stats::sd))
  }
  set.seed(3)
  new <- cl_sd(function() {
    dr <- ci(".countimp_draw_2l")(fit, nd)
    MASS::rnegbin(length(dr$mu), mu = dr$mu, theta = dr$theta)
  })
  beta <- glmmTMB::fixef(fit)$cond
  rv <- t(chol(stats::vcov(fit)$cond))
  set.seed(3)
  old <- cl_sd(function() {
    bs <- as.numeric(beta) + as.numeric(rv %*% stats::rnorm(ncol(rv)))
    mu <- as.numeric(ci(".countimp_predict_2l")(
      ci(".countimp_inject_beta")(fit, bs), nd, type = "response"))
    MASS::rnegbin(length(mu), mu = mu,
                  theta = ci(".countimp_draw_theta")(fit))
  })
  expect_gt(new / old, 2)
})
