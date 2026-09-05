test_that(".cm_detect_conflicts returns nothing with fewer than two attached packages", {
  expect_identical(.cm_detect_conflicts("notattached_pkg_xyz"), character(0))
  expect_identical(.cm_detect_conflicts(c("stats", "notattached_pkg_xyz")), character(0))
})

test_that(".cm_detect_conflicts reports names exported by several attached packages", {
  skip_if_not("package:dplyr" %in% search(), "dplyr is not attached")
  out <- .cm_detect_conflicts(c("stats", "dplyr"))
  expect_true(any(grepl("^filter\\(\\): dplyr masks stats$", out)))
  expect_true(any(grepl("^lag\\(\\): dplyr masks stats$", out)))
})

test_that(".cm_detect_conflicts hides the routine data.table/tidyverse overlaps", {
  skip_if_not("package:dplyr" %in% search(), "dplyr is not attached")
  skip_if_not("package:data.table" %in% search(), "data.table is not attached")
  out <- .cm_detect_conflicts(c("data.table", "dplyr"))
  expect_false(any(grepl("^(first|last|between)\\(\\)", out)))
})
