# countimp 3.0.1

**Behaviour change: the imputations of every `2l.*.boot` method are different
from 3.0.0.** The cluster bootstrap treated a cluster the resample happened to
leave out as an unseen cluster and drew its effect from the prior instead of
its own conditional posterior -- about a third of the imputed rows, from a
distribution roughly five times too wide. The Bayesian variants are unchanged,
bit for bit. Re-run any study whose bootstrap arms were computed with 3.0.0;
the coverage of the dispersion parameter rose from 31.5 % to 96.5 % in ours.

* **`2l.*` methods work again through `method` + `predictorMatrix`.** The level
  check indexed the method vector by name, so the unnamed form that mice's own
  documentation uses -- `method = c("2l.nb2", "", "", "")` -- aborted with
  "subscript out of bounds" before anything ran. Named partial vectors,
  `method = c(y = "2l.nb2")`, silently selected the wrong variables. Both forms
  are now normalised onto the columns of the predictorMatrix.

* **A grouping variable may be a factor on the formula route.** `(1 | g)` is
  lme4's and glmmTMB's syntax, where a factor is the normal case. The classic
  route still requires an integer -- there a factor could be expanded into
  indicator columns -- but now names the variable and gives the fix.

* **A `2lonly.*` method aimed at too fine a level is refused.** Grouping a
  school-level variable by class passes every check the aggregation itself can
  make and then draws one value per class, leaving a school with several
  different values. It is an error now, not silent damage.

* Clearer messages where the previous ones misled: an underdispersed negative
  binomial no longer blames collinear predictors; a singular fixed-effect
  covariance names the situation rather than reporting a LAPACK minor.

# countimp 3.0.0

## New

* **A vignette**: `vignette("countimp")` -- choosing a model family, the
  formula interface, why a log-link model with count predictors extrapolates,
  and how to check imputations.

* **Bounded, censored and underdispersed counts.** `bp`/`bnb` for a variable
  with a known range (`bounds = c(0, 10)`), `cp`/`cnb` for right-censored
  counts (`censor = 10`), `cmp` for Conway-Maxwell-Poisson, which is the only
  family here that also covers UNDERdispersion. All single-level; the scale
  argument is required and is not guessed from the data.

* **Zero-truncated counts**: `ztp`/`ztnb` for variables that cannot be zero --
  days in hospital among those admitted, offences among offenders.

* **Cluster-level variables**: `2lonly.pmm`, `2lonly.norm`, `2lonly.pois`,
  `2lonly.nb` aggregate to one row per cluster, impute there, and hand the
  value back to every row of it. Where several grouping columns are coded
  `-2`, the outermost one the variable is constant within is used, so three
  levels need no separate method.

* **Two-level COM-Poisson**: `2l.cmp` and `2l.cmp.boot`.

## Behaviour changes

* **Imputing a cluster-level variable row by row is refused.** A variable that
  is constant within its cluster and imputed per row yields completed data the
  design cannot contain. Use a `2lonly.*` method, or
  `options(countimp.check.levels = FALSE)` to override.

* **An observed zero under a zero-truncated method is an error**, not a
  warning: the value lies outside the model's support.

* **Giving a scale argument both through the family and as an argument is an
  error**, not a precedence rule -- the two can differ and the imputations
  would not show it.

* **A named `method` vector is matched by name**, not by position, and naming
  only some variables no longer leaves the others unimputed.

* The `.boot` variants require nested grouping factors; a crossed design is
  refused rather than silently resampled.

## Fixes

* A model meant as three-level ran as two-level, without a warning.
* The dispersion is no longer drawn where the data do not determine it; a
  draw far above everything observed is reported.
* Wrong-family messages name the actual problem -- `family = poisson()` where
  a countimp family is expected used to fail deep inside the engine.
* `MASS` was listed in both Imports and Suggests.

# countimp 2.1.0 - 2.6.0 (internal; never released)

The last public release was 2.0.7. Everything between it and version 3 was a
development step that nobody outside the project ran, so this section says
what changed between 2.0.7 and 3:

* **Automatic outlier handling is off by default (2.6.0).** `EV = TRUE` now
  raises a deprecation warning.

* **The dispersion is no longer drawn where the data do not determine it
  (2.6.0).** A negative binomial fitted to data without overdispersion drives
  theta towards infinity and its standard error with it; drawing from that
  produced imputations orders of magnitude above anything observed. Detected
  by max SE(beta) > 5 or theta < 0.25, then refitted with theta fixed.

* **'aster' is no longer required (2.3.0).** The zero-truncated draws are
  implemented in base R, verified against aster 1.3.7 across 17 parameter
  combinations. Installing countimp no longer pulls in 3.1 MB of compiled
  code for two functions. The base R draws are about ten times slower, which
  does not matter at the sizes involved.

* **A diagnostics layer (2.1.0).** `countimp_diagnostics()` records what each
  draw did; the checks FLAG, they do not silently repair. Introduced after
  deprecating `EV = TRUE` removed a safeguard without putting anything in its
  place.

