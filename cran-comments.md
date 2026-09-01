## Submission

New submission.

`countimp` has been developed and distributed via GitHub
(<https://github.com/kkleinke/countimp>) since 2013 and is used in the applied
missing-data literature; the version number therefore continues the existing
GitHub series rather than starting at 0.1.0. Version 3.0.0 is the first version
prepared for CRAN: the imputation engine is now self-contained, so the package
no longer depends on `mice`, and imputation models can be specified through a
formula interface.

## Test environments

* local macOS (arm64), R 4.x -- `R CMD check --as-cran`
* R-hub and win-builder: to be run before submission

## R CMD check results

0 errors | 0 warnings | 3 notes

* `New submission` -- see above.
* `checking for future file timestamps ... unable to verify current time`
  -- the check machine had no network access to the time server; not a
  package issue.
* `Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being
  installed` -- pandoc is absent on the local check machine.

## Notes for the reviewer

* `mice` is in `Suggests`, not `Imports`. The package's own engine performs the
  imputations; the only `mice::` call in the code is a fallback for `mids`
  objects the package cannot read, guarded by `.countimp_need_mice()`, which
  wraps `requireNamespace("mice", quietly = TRUE)` and stops with an
  instructive message. Verified with `.Library` redirected to a `mice`-free
  library: `countimp()`, `countimp_complete()`, `with()` and `miinference()`
  all complete. Tests that need `mice` call
  `skip_if_not_installed("mice")`.
* The 63 exported `mice.impute.*` functions are the documented method names of
  the classic interface and are retained for backwards compatibility. Since
  3.0.0 they need not be chosen by hand -- the formula interface selects them.
* Examples that fit multilevel count models are wrapped in `\donttest{}` where
  they exceed the 5-second budget.
