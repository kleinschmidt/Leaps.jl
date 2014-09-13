# Leaps: Regression subset selection

This package performs (or will eventually) subset selection for regression predictors, using forward, backard, and exhaustive search.  Currently exhaustive search is implemented, using the `xhaust!(::RegSubsets)` command.

It is based on the `leaps` package for `R` by Thomas Lumley (<tlumley@u.washington.edu>), which is in turn based on Fortran code wraps the Fortran code for subset selection by by Alan Miller of CSIRO Division of Mathematics & Statistics (<Alan.Miller@vic.cmis.csiro.au>).

