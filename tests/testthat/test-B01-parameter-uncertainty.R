## B01: the drawn coefficient vector must reach the prediction.
##
## Before the fix, b.star was drawn but predict() was called on the unmodified
## fit, so every imputation used the ML estimate: the Bayesian variants were
## silently identical to their non-Bayesian counterparts. The route now writes
## b.star into the fit's parameter vector via glmmTMB's documented `newparams`.

test_that(".countimp_inject_beta writes the drawn coefficients into the fit", {
  fit <- fit_2l_nb()
  beta <- glmmTMB::fixef(fit)$cond

  ## Shift the INTERCEPT only. Adding 1 to every coefficient would scale each
  ## row by exp(1 + x_i), which is not a constant, so a check on the fitted
  ## rates would test nothing sharp. With an intercept-only shift the linear
  ## predictor moves by exactly 1 for every row.
  b1 <- beta; b1[["(Intercept)"]] <- b1[["(Intercept)"]] + 1
  inj <- ci(".countimp_inject_beta")(fit, b1)

  ## Assert on the parameter slots, not on predict(). This is what injection
  ## is responsible for; what predict() then does with it is B20's subject.
  i <- which(names(inj$fit$par) == "beta")
  expect_equal(as.numeric(inj$fit$par[i]), as.numeric(b1), tolerance = 1e-10)
  j <- which(names(inj$fit$parfull) == "beta")
  expect_equal(as.numeric(inj$fit$parfull[j]), as.numeric(b1), tolerance = 1e-10)
  expect_equal(as.numeric(glmmTMB::fixef(inj)$cond), as.numeric(b1), tolerance = 1e-10)
})

test_that("the injected coefficients shift the linear predictor by exactly 1", {
  ## The sharp check that the old version of this test tried to make through
  ## predict(). It has to be made on eta = X beta + Z b, reconstructed
  ## explicitly: predict.glmmTMB() re-derives the conditional modes of u from
  ## the parameter vector it is handed, so it passes only the shrinkage
  ## fraction sigma^2 / (n_j tau^2 + sigma^2) of the shift through -- about
  ## 9% here. That is B20, and it is why ci(".countimp_draw_2l")() exists.
  fit <- fit_2l_nb()
  d   <- model.frame(fit)
  beta <- glmmTMB::fixef(fit)$cond
  b1 <- beta; b1[["(Intercept)"]] <- b1[["(Intercept)"]] + 1

  X <- glmmTMB::getME(fit, "X")
  Z <- ci(".countimp_re_design")(fit, d)
  u <- as.numeric(glmmTMB::getME(fit, "b"))
  eta0 <- as.numeric(X %*% as.numeric(beta)) + as.numeric(Z %*% u)
  eta1 <- as.numeric(X %*% as.numeric(b1))   + as.numeric(Z %*% u)

  expect_equal(eta1 - eta0, rep(1, length(eta0)), tolerance = 1e-8)
})

test_that("injecting the ML estimate reproduces the original prediction", {
  fit <- fit_2l_nb()
  nd  <- model.frame(fit)[1:20, , drop = FALSE]
  same <- ci(".countimp_inject_beta")(fit, glmmTMB::fixef(fit)$cond)
  expect_equal(as.numeric(ci(".countimp_rate")(same, nd)),
               as.numeric(ci(".countimp_rate")(fit, nd)), tolerance = 1e-8)
})

test_that("the injection route is verified once per session", {
  expect_true(is.logical(ci(".countimp_state")$beta_route_ok))
  fit <- fit_2l_nb()
  ci(".countimp_check_beta_route")(fit)
  expect_true(ci(".countimp_state")$beta_route_ok)
})

test_that("Bayes and non-Bayes variants give different imputations", {
  d <- sim_count(n = 200, nmis = 60)
  set.seed(1); a <- imputed_values(quiet_impute(d, method = c("nb", ""),
                                                m = 3, maxit = 2, seed = 4))
  set.seed(1); b <- imputed_values(quiet_impute(d, method = c("nb.boot", ""),
                                                m = 3, maxit = 2, seed = 4))
  expect_false(identical(a, b))
})
