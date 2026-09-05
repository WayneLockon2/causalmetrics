test_that(".cm_as_dt coerces data.frame and passes data.table through", {
  df <- data.frame(a = 1:3, b = letters[1:3])
  out <- .cm_as_dt(df)
  expect_s3_class(out, "data.table")
  expect_equal(nrow(out), 3L)

  dt <- data.table::data.table(a = 1:3)
  expect_identical(.cm_as_dt(dt), dt)
})

test_that(".cm_as_dt rejects non-data.frame input", {
  expect_error(.cm_as_dt(1:3), "data.frame or data.table")
  expect_error(.cm_as_dt(list(a = 1)), "data.frame or data.table")
})

test_that(".cm_assert_cols catches missing columns", {
  dt <- data.table::data.table(a = 1, b = 2)
  expect_null(.cm_assert_cols(dt, c("a", "b")))
  expect_error(.cm_assert_cols(dt, c("a", "zzz")), "Missing columns: zzz")
  expect_error(.cm_assert_cols(data.frame(a = 1), "a"))
})

test_that(".cm_check_package gives an informative error", {
  expect_error(
    .cm_check_package("notARealPackage12345"),
    "Package 'notARealPackage12345' is required"
  )
  expect_silent(.cm_check_package("stats"))
})
