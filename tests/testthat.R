# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

if (requireNamespace("testthat", quietly = TRUE)) {
  testthat::test_check("rbims")
}


library(testthat)
library(rbims)

test_check("rbims")
