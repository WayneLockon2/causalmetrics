test_that("table_task returns a booktabs LaTeX kable", {
  skip_if_not_installed("kableExtra")
  dt <- data.frame(group = c("A", "B"), mean = c(1.234, 5.678))
  tbl <- table_task(dt, digits = 2)
  expect_s3_class(tbl, "knitr_kable")
  expect_equal(attr(tbl, "format"), "latex")
  tex <- as.character(tbl)
  expect_match(tex, "\\\\toprule")
  expect_match(tex, "1.23", fixed = TRUE)
})

test_that("wrap_latex_table builds a complete table environment", {
  skip_if_not_installed("kableExtra")
  tbl <- table_task(data.frame(x = 1:2, y = c(0.5, 1.5)))

  out <- wrap_latex_table(tbl, caption = "Cap", notes = "Some notes.",
                          label = "tab:x", width = 0.8)
  expect_s3_class(out, "cm_latex_block")
  tex <- as.character(out)
  expect_match(tex, "^\\\\begin\\{table\\}\\[H\\]")
  expect_match(tex, "\\\\caption\\{Cap\\}\\\\label\\{tab:x\\}")
  expect_match(tex, "\\\\resizebox\\{0.8\\\\textwidth\\}\\{!\\}")
  expect_match(tex, "\\\\textit\\{Notes\\}: Some notes\\.")
  expect_match(tex, "\\\\end\\{table\\}$")
  # Exactly one table environment: the kable's own wrapper is stripped.
  expect_equal(lengths(regmatches(tex, gregexpr("\\\\begin\\{table\\}", tex))), 1L)

  out2 <- wrap_latex_table(tbl, caption = "Cap")
  expect_false(grepl("resizebox", as.character(out2)))
  expect_false(grepl("\\\\label", as.character(out2)))
  expect_false(grepl("Notes", as.character(out2)))

  out3 <- wrap_latex_table(tbl, width = 0.5, height = 0.3)
  expect_match(as.character(out3), "\\\\resizebox\\{0.5\\\\textwidth\\}\\{0.3\\\\textheight\\}")

  expect_error(wrap_latex_table(tbl, width = 2), "width")
  expect_error(wrap_latex_table(tbl, caption = 1), "caption")
})

test_that("wrap_latex_table can write the tabular to a file and input it", {
  skip_if_not_installed("kableExtra")
  path <- file.path(tempfile("tbl"), "inner.tex")
  out <- wrap_latex_table(table_task(data.frame(a = 1)), caption = "C", path = path)
  expect_true(file.exists(path))
  expect_true(grepl(sprintf("\\input{%s}", path), as.character(out), fixed = TRUE))
  expect_match(paste(readLines(path), collapse = "\n"), "tabular")
})

test_that(".cm_strip_outer_table_env removes the outer environment, centering and caption", {
  tex <- paste(
    "\\begin{table}[H]", "\\centering", "\\caption{Old}\\label{tab:old}",
    "\\begin{tabular}{l}", "x\\\\", "\\end{tabular}", "\\end{table}",
    sep = "\n"
  )
  inner <- .cm_strip_outer_table_env(tex)
  expect_false(grepl("begin\\{table\\}", inner))
  expect_false(grepl("centering", inner))
  expect_false(grepl("caption", inner))
  expect_match(inner, "^\\\\begin\\{tabular\\}")

  bare <- "\\begin{tabular}{l}x\\end{tabular}"
  expect_identical(.cm_strip_outer_table_env(bare), bare)
})

test_that("kable_notes appends a notes line and keeps the kable class", {
  skip_if_not_installed("kableExtra")
  tbl <- table_task(data.frame(a = 1))
  out <- kable_notes(tbl, "Note text.")
  expect_s3_class(out, "knitr_kable")
  expect_match(as.character(out), "\\\\textit\\{Notes\\}: Note text\\.")
  expect_error(kable_notes(tbl, notes = 1), "string")
})

test_that("cm_latex_block prints raw LaTeX and knits as-is", {
  blk <- structure("\\begin{table}\\end{table}",
                   class = c("cm_latex_block", "character"), format = "latex")
  expect_output(print(blk), "\\\\begin\\{table\\}")
  res <- NULL
  capture.output(res <- withVisible(print(blk)))
  expect_false(res$visible)
  expect_identical(res$value, blk)
  expect_s3_class(knitr::knit_print(blk), "knit_asis")
})
