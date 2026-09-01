# Extracted from test-B04-underdispersion.R:32

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "countimp", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
st <- ci(".countimp_state")
st$underdisp_warned <- FALSE
d <- sim_underdispersed()
w <- character(0)
withCallingHandlers(
    countimp(d, method = c("quasipoisson", ""), m = 3, maxit = 3, seed = 5,
             printFlag = FALSE),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
    })
expect_equal(sum(grepl("underdispers", w)), 1L)
