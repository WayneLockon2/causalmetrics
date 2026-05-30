# R/plot-task.R

#' Standardized plot template
#'
#' Wraps `ggplot()` with the project's standard theme, color scales, and
#' legend position. Use `wrap_latex_figure()` after composing the plot to
#' produce a LaTeX figure block for inclusion in a paper.
#'
#' @param ... Arguments passed to `ggplot()`. Typically `data` and `mapping`.
#'   Add layers with `+` exactly as you would for `ggplot()`.
#'
#' @return A ggplot object with the standard causalmetrics styling applied.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- plot_task(mtcars, aes(wt, mpg, color = factor(cyl))) +
#'   geom_point() +
#'   labs(x = "Weight", y = "MPG")
#'
#' # For paper-ready LaTeX output:
#' tex <- wrap_latex_figure(
#'   p,
#'   caption = "Fuel efficiency vs. weight",
#'   notes   = "Each point is one car.",
#'   label   = "fig:mpg_weight",
#'   path    = "figures/mpg_weight.pdf"
#' )
#' }
#'
#' @export
plot_task <- function(...) {

    base <- ggplot2::ggplot(...) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::scale_linetype_manual(
            values = c("solid", "longdash", "dotted",
                       "dashed", "dotdash", "twodash")
        )

    if (requireNamespace("ggthemes", quietly = TRUE)) {
        base <- base +
            ggthemes::scale_colour_stata("s2color") +
            ggthemes::scale_fill_stata("s2color")
    }

    base
}


#' Save a ggplot and return a LaTeX figure block
#'
#' Takes a composed ggplot, saves it to PDF, and returns a
#' `\\begin{figure}` block referencing the saved file. Use after
#' `plot_task()` once you've added all geoms and labels.
#'
#' @param plot A ggplot object.
#' @param caption Character. Figure caption.
#' @param notes Character. Italicized notes line below the figure.
#' @param label Character. LaTeX label, e.g., `"fig:event_study"`.
#' @param width Numeric. Figure width as a fraction of `\\textwidth`.
#'   Default 0.8.
#' @param path Character or NULL. File path for the saved PDF. If `NULL`
#'   (default), a tempfile in `tempdir()` is used. Tempfile output is fine
#'   for one-shot rendering in the same R session, but for documents you
#'   plan to re-knit later, pass an explicit path like
#'   `"figures/event_study.pdf"`.
#' @param fig_width Numeric. PDF width in inches. Default 7.
#' @param fig_height Numeric. PDF height in inches. Default 4.5.
#'
#' @return Character string of LaTeX code (with class `cm_latex_block`).
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- plot_task(mtcars, aes(wt, mpg)) + geom_point()
#'
#' # Quick render to tempfile
#' wrap_latex_figure(p, caption = "Quick look", label = "fig:quick")
#'
#' # Persistent path for a paper
#' wrap_latex_figure(
#'   p,
#'   caption = "Final figure",
#'   label   = "fig:final",
#'   path    = "figures/final.pdf"
#' )
#' }
#'
#' @export
wrap_latex_figure <- function(plot,
                              caption    = NULL,
                              notes      = NULL,
                              label      = NULL,
                              width      = 0.8,
                              path       = NULL,
                              fig_width  = 7,
                              fig_height = 4.5) {

    checkmate::assert_class(plot, "ggplot")
    checkmate::assert_string(caption, null.ok = TRUE)
    checkmate::assert_string(notes,   null.ok = TRUE)
    checkmate::assert_string(label,   null.ok = TRUE)
    checkmate::assert_number(width,  lower = 0, upper = 1)
    checkmate::assert_number(fig_width,  lower = 0.1)
    checkmate::assert_number(fig_height, lower = 0.1)
    checkmate::assert_string(path, null.ok = TRUE)

    # Default to a tempfile under tempdir(). R cleans up tempdir on session
    # exit, so the PDF persists long enough to compile the document but does
    # not pollute the user's working directory.
    if (is.null(path)) {
        path <- tempfile(pattern = "cm_fig_", tmpdir = tempdir(),
                         fileext = ".pdf")
    }

    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

    ggplot2::ggsave(
        filename = path,
        plot     = plot,
        width    = fig_width,
        height   = fig_height,
        device   = "pdf"
    )

    cap_line <- if (!is.null(caption)) {
        if (!is.null(label)) {
            sprintf("\\caption{%s}\\label{%s}", caption, label)
        } else {
            sprintf("\\caption{%s}", caption)
        }
    } else ""

    notes_line <- if (!is.null(notes)) {
        sprintf(
            "\\par\\vspace{0.5em}\\noindent{\\footnotesize \\textit{Notes}: %s}",
            notes
        )
    } else ""

    tex <- paste(
        "\\begin{figure}[H]",
        "\\centering",
        cap_line,
        sprintf("\\includegraphics[width=%g\\textwidth]{%s}", width, path),
        notes_line,
        "\\end{figure}",
        sep = "\n"
    )

    structure(
        tex,
        class  = c("cm_latex_block", "character"),
        format = "latex"
    )
}
