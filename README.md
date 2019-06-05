countimp: Multiple Imputation of incomplete count data
======================================================

Installation
------------

The latest version can be installed from GitHub using function `install_git` from package `devtools`:

``` r
devtools::install_git(url = "https://github.com/kkleinke/countimp", 
                      branch = "master")
```

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
-------------------------------------------------------------------------

References
----------
Kleinke, K., & Reinecke, J. (2013a). Multiple imputation of incomplete zero-inflated count data. Statistica Neerlandica, 67(3), 311–336. doi: 10.1111/stan.12009
http://onlinelibrary.wiley.com/doi/10.1111/stan.12009/abstract

Kleinke, K., & Reinecke, J. (2015a). Multiple imputation of multilevel count data. In U. Engel, B. Jann, P. Lynn, A. Scherpenzeel, and P. Sturgis (Eds.), Improving Survey Methods: Lessons from Recent Research (pp. 381–396). New York: Routledge, Taylor & Francis Group.
http://www.psypress.com/books/details/9780415817622/

Kleinke, K., & Reinecke, J. (2015b). Multiple imputation of overdispersed multilevel count data. In: Uwe Engel (Ed.), Survey Measurements. Techniques, Data Quality and Sources of Error (pp. 209–226). Frankfurt a. M.: Campus/The University of Chicago Press.
http://press.uchicago.edu/ucp/books/book/distributed/S/bo22196267.html

Further information
-------------------
https://www.kkleinke.com/software
