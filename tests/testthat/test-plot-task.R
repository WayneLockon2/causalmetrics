test_that("plot_task returns a themed ggplot that builds", {
  p <- plot_task(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
  expect_s3_class(p, "ggplot")
  expect_silent(built <- ggplot2::ggplot_build(p))
  expect_equal(nrow(built$data[[1]]), nrow(mtcars))
})

test_that("wrap_latex_figure saves a PDF and returns a figure block", {
  p <- plot_task(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
  path <- file.path(tempfile("figdir"), "plot.pdf")

  out <- wrap_latex_figure(p, caption = "Cap", notes = "N.", label = "fig:x",
                           width = 0.6, path = path)
  expect_s3_class(out, "cm_latex_block")
  expect_true(file.exists(path))
  expect_gt(file.size(path), 0)

  tex <- as.character(out)
  expect_match(tex, "^\\\\begin\\{figure\\}\\[H\\]")
  expect_match(tex, "\\\\caption\\{Cap\\}\\\\label\\{fig:x\\}")
  expect_match(tex, "\\\\includegraphics\\[width=0.6\\\\textwidth\\]")
  expect_true(grepl(path, tex, fixed = TRUE))
  expect_match(tex, "\\\\textit\\{Notes\\}: N\\.")
  expect_match(tex, "\\\\end\\{figure\\}$")

  out2 <- wrap_latex_figure(p)
  expect_match(as.character(out2), "cm_fig_")
  expect_false(grepl("caption", as.character(out2)))

  expect_error(wrap_latex_figure("not a plot"), "ggplot")
  expect_error(wrap_latex_figure(p, width = 1.5), "width")
})
