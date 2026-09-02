countimp: Multiple Imputation of Incomplete Count Data
======================================================

Further information
-------------------
[https://kristian-kleinke.de/appliedmi.html#r-package-countimp---multiple-imputation-of-incomplete-count-data](https://kristian-kleinke.de/appliedmi.html#r-package-countimp---multiple-imputation-of-incomplete-count-data)

User's Guide
-------------------

[https://kristian-kleinke.de/countimp](https://kristian-kleinke.de/countimp/)

The guide and the manual still describe version 2.x and are being revised. For
version 3, the help pages and `vignette("countimp")` are the current reference.

Installation
------------

`countimp` is not on CRAN; install it from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("kkleinke/countimp")
```

The package brings its own imputation engine and needs no compiler. `mice` is
*optional* -- install it only if you also need non-count imputation methods or
`mice`'s pooling functions.

If loading `countimp` warns that `glmmTMB was built with TMB package version
...`, the installed `glmmTMB` binary was compiled against an older `TMB` than
the one you now have. It comes from `glmmTMB`, not from `countimp`, and
rebuilding it from source resolves it:

``` r
install.packages("glmmTMB", type = "source")
```

Specifying imputation models
----------------------------

Since version 3.0.0 there are two ways to specify the imputation models. The
formula interface is the recommended one:

``` r
library(countimp)
data(crim4w)

## one formula per incomplete variable; one family covers all of them
imp <- countimp(data     = crim4w,
                formulas = list(ACRIM ~ FEMALE + RE + GY + HA,
                                BCRIM ~ FEMALE + RE + GY + HA,
                                CCRIM ~ FEMALE + RE + GY + HA,
                                DCRIM ~ FEMALE + RE + GY + HA),
                family   = nb(),
                m = 5, maxit = 5, seed = 1)
```

The same four models through the classic interface -- a predictor matrix and
one method name per variable:

``` r
pred <- countimp(crim4w, maxit = 0)$predictorMatrix
pred[, ] <- 0
pred[c("ACRIM", "BCRIM", "CCRIM", "DCRIM"),
     c("FEMALE", "RE", "GY", "HA")] <- 1

imp <- countimp(data            = crim4w,
                method          = c(ACRIM = "nb", BCRIM = "nb",
                                    CCRIM = "nb", DCRIM = "nb"),
                predictorMatrix = pred,
                m = 5, maxit = 5, seed = 1)
```

Both calls give the same imputations, bit for bit, from the same seed. The
difference is where the model lives: the formula carries the predictors of one
variable next to that variable, while the matrix states them for all variables
at once -- here by zeroing every cell and switching sixteen of them back on.
The default fills the matrix with every other column, which for these data
means the four counts predict each other and `id`, a person identifier, enters
as a predictor as well.

Predictor roles and the random-effects structure are read from the formula, so
no `predictorMatrix` and no hand-picked method name are needed -- 70 of the 85
`mice.impute.*` variants are selected for you. To see how a specification
is translated before imputing:

``` r
data(crim4l)   # long format, ID identifies the person
countimp_spec(data     = crim4l,
              formulas = list(DELINQ = DELINQ ~ FEMALE + TIME + (1 | ID)),
              family   = list(DELINQ = hurdle_nb()),
              draw     = "boot")$method["DELINQ"]
#>        DELINQ
#> "2l.hnb.boot"
```

The basic family constructors are `poisson_count()`, `nb()`, `zi_poisson()`,
`zi_nb()`, `hurdle_poisson()`, `hurdle_nb()` and `compois()`; `draw` is
`"bayes"` (default) or `"boot"`. Two-part families take a separate `zero`
formula for the zero/hurdle part.

**Counts with a restricted support** have their own families, and the
restriction is part of the family rather than of the formula:

``` r
countimp(d,
         formulas = list(visits ~ x1, days ~ x1, spells ~ x1),
         family   = list(visits = censored_poisson(censor = 10),  # top-coded
                         days   = bounded_poisson(bounds = c(0, 5)),
                         spells = zerotrunc_poisson()),           # never zero
         m = 5)
```

The scale then stands next to the variable it describes. `bounds` and `censor`
are mandatory in these families; giving the same scale both through the family
and as a separate argument is an error, not a precedence rule. They exist for
single-level data only -- there is no `2l.bp` or `2l.cp`, and a grouping term
in the formula is refused rather than compiled into a method that does not
exist.

The classic interface remains unchanged and is the way to supply a scale
separately:

``` r
countimp(d, method = c(visits = "cp", x1 = ""), m = 5, censor = 10)
```

**Three levels.** A formula may carry more than one grouping term, which is how
a three-level model is written -- pupils in classes in schools:

``` r
countimp(d, formulas = list(y ~ x1 + (1 | school) + (1 | class)),
         family = poisson_count(), m = 5)
```

No new method name is involved: the levels travel in the type codes, so
`2l.poisson` and its relatives fit two, three or more levels. The first
grouping term carries the random slopes; further terms enter as random
intercepts. The `.boot` variants are two-level only for now -- a nested design
needs hierarchical resampling, and they say so rather than resampling the outer
level alone.

The two interfaces are mutually exclusive: `formulas` together with `method` or
`predictorMatrix` is an error.

The classic interface -- a `predictorMatrix`, type codes and an explicit method
name per variable -- continues to work unchanged.

Main functions
--------------

| Function Name                | model / family                         |
|------------------------------|----------------------------------------|
| `mice.impute.poisson()`      | Poisson                                |
| `mice.impute.quasipoisson()` | Quasi-Poisson                          |
| `mice.impute.nb()`           | negative binomial                      |
| `mice.impute.zip()`          | zero-inflated Poisson                  |
| `mice.impute.zinb()`         | zero-inflated negative binomial        |
| `mice.impute.hp()`           | hurdle Poisson                         |
| `mice.impute.hnb()`          | hurdle negative binomial               | 
| `mice.impute.2l.poisson()`   | two-level Poisson                      |
| `mice.impute.2l.nb()`        | two-level negative binomial            |
| `mice.impute.2l.zip()`       | two-level zero-inflated Poisson        |
| `mice.impute.2l.zinb()`      | two-level zero-inflated NB             |
| `mice.impute.2l.hp()`        | two-level hurdle Poisson               |
| `mice.impute.2l.hnb()`       | two-level hurdle negative binomial     |
| `mice.impute.ztp()`          | zero-truncated Poisson                 |
| `mice.impute.ztnb()`         | zero-truncated negative binomial       |
| `mice.impute.bp()`           | bounded Poisson (`bounds`)             |
| `mice.impute.bnb()`          | bounded negative binomial (`bounds`)   |
| `mice.impute.cp()`           | right-censored Poisson (`censor`)      |
| `mice.impute.cnb()`          | right-censored negative binomial (`censor`) |
| `mice.impute.cmp()`          | Conway-Maxwell-Poisson (underdispersion) |
| `mice.impute.2l.cmp()`       | two-level Conway-Maxwell-Poisson       |
| `mice.impute.2l.pmm()`       | two-level predictive mean matching     |
| `mice.impute.2lonly.pois()`  | a count variable that is constant within its cluster |
| `mice.impute.2lonly.nb()`    | the same, negative binomial            |
| `mice.impute.2lonly.pmm()`   | the same, predictive mean matching     |
| `mice.impute.2lonly.norm()`  | the same, normal                       |
-------------------------------------------------------------------------

The `2lonly` methods are for variables that live at the cluster level -- class
size, ward type, school budget. They aggregate to one row per cluster, impute
there and hand the value back to every row of it, so a class cannot end up
with two different sizes. Code the grouping column as `-2` in the
`predictorMatrix`. With more than one such column -- three-level data -- the
outermost level the variable is constant within is used: a school-level
variable is constant within classes too, and imputing it per class would give
classes of the same school different values.

Every family above has a `.boot` variant that draws its parameters from a
bootstrap instead of the normal approximation -- a cluster bootstrap for the
two-level ones; the exceptions are `2l.pmm` and the four `2lonly` methods.

A `.noint` variant holds the intercept fixed across clusters. `2l.poisson`,
`2l.nb` and `2l.cmp` spell it `.noint`; the two-part models `2l.zip`, `2l.zinb`,
`2l.hp` and `2l.hnb` say which part is meant -- `.noint.count`, `.noint.zero`
or `.noint.both`. The two can be combined: `mice.impute.2l.zip.noint.both.boot`.

`countimp_fit_diag()` compares the candidate families against an observed
count variable and recommends one.

References
----------
van Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate Imputation by Chained Equations in R. Journal of Statistical Software, 45(3), 1–67. doi: 10.18637/jss.v045.i03

Kleinke, K., & Reinecke, J. (2013a). Multiple imputation of incomplete zero-inflated count data. Statistica Neerlandica, 67(3), 311–336. doi: 10.1111/stan.12009

Kleinke, K., & Reinecke, J. (2015a). Multiple imputation of multilevel count data. In U. Engel, B. Jann, P. Lynn, A. Scherpenzeel, and P. Sturgis (Eds.), Improving Survey Methods: Lessons from Recent Research (pp. 381–396). New York: Routledge, Taylor & Francis Group.
http://www.psypress.com/books/details/9780415817622/

Kleinke, K., & Reinecke, J. (2015b). Multiple imputation of overdispersed multilevel count data. In: Uwe Engel (Ed.), Survey Measurements. Techniques, Data Quality and Sources of Error (pp. 209–226). Frankfurt a. M.: Campus/The University of Chicago Press.
http://press.uchicago.edu/ucp/books/book/distributed/S/bo22196267.html

Kleinke, K., Reinecke, J., Salfrán, D., & Spiess, M. (2020). Applied Multiple Imputation: Advantages, Pitfalls, New Developments and Applications in R. Cham: Springer. doi: 10.1007/978-3-030-38164-6
https://doi.org/10.1007/978-3-030-38164-6

