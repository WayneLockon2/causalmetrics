# R/table-task.R

#' Standardized table template
#'
#' Wraps `knitr::kable()` with `booktabs` formatting and `kableExtra::kable_styling()`
#' defaults. Returns a kable object suitable for inline RMarkdown rendering.
#' Use `wrap_latex_table()` to wrap the result in a full `\\begin{table}` block
#' with caption, notes, and resizebox for paper output.
#'
#' @param x A data.frame, data.table, or matrix.
#' @param ... Additional arguments passed to `knitr::kable()`. E.g., `digits`,
#'   `col.names`, `align`, `escape`.
#'
#' @return A kable/kable_styling object.
#'
#' @examples
#' \dontrun{
#' library(data.table)
#' dt <- data.table(group = c("A", "B"), mean = c(1.2, 2.4), n = c(100, 150))
#'
#' # Inline RMarkdown rendering:
#' table_task(dt, digits = 2)
#'
#' # For paper-ready LaTeX:
#' tbl <- table_task(dt, digits = 2)
#' tex <- wrap_latex_table(
#'   tbl,
#'   caption = "Group summary statistics",
#'   notes   = "Means and observation counts by group.",
#'   label   = "tab:summary",
#'   width   = 0.6
#' )
#' writeLines(tex, "tables/summary.tex")
#' }
#'
#' @export
table_task <- function(x, ...) {
    .cm_check_package("knitr")
    .cm_check_package("kableExtra")

    tbl <- knitr::kable(
        x,
        booktabs = TRUE,
        format   = "latex",
        linesep  = "",
        ...
    )

    kableExtra::kable_styling(
        tbl,
        position      = "center",
        latex_options = "HOLD_position"
    )
}

#' Wrap a kable in a LaTeX table environment
#'
#' Takes a finished kable and produces a complete `\\begin{table}` block with
#' caption, optional resizebox, and a notes line styled as
#' `{\\footnotesize \\textit{Notes}: ...}`. Use after `table_task()` once the
#' table content is finalized.
#'
#' @param tbl A kable or kable_styling object, ideally from `table_task()`.
#' @param caption Character. Table caption. Default `NULL`.
#' @param notes Character. Italicized notes line below the table.
#' @param label Character. LaTeX label, e.g., `"tab:summary"`.
#' @param width Numeric or NULL. Table width as a fraction of `\\textwidth`.
#'   If `NULL`, no resizebox is applied and the table renders at its natural
#'   width. Default `NULL`.
#' @param height Numeric or NULL. Table height as a fraction of
#'   `\\textheight`. Used only if `width` is provided. Default `NULL`
#'   (preserves aspect from width).
#' @param path Character or NULL. If provided, the inner tabular is written
#'   to this `.tex` file and the LaTeX block uses `\\input{path}`. If `NULL`,
#'   the tabular is embedded inline. Default `NULL`.
#'
#' @return Character string of LaTeX code.
#'
#' @examples
#' \dontrun{
#' dt <- data.frame(x = 1:3, y = c(2.1, 3.4, 4.5))
#' tbl <- table_task(dt, digits = 2)
#'
#' # Simple wrap, no resizebox:
#' wrap_latex_table(tbl, caption = "Example", label = "tab:ex")
#'
#' # With resizebox and notes:
#' wrap_latex_table(
#'   tbl,
#'   caption = "Example",
#'   notes   = "Some explanation here.",
#'   label   = "tab:ex",
#'   width   = 0.8
#' )
#' }
#'
#' @export
wrap_latex_table <- function(tbl,
                             caption = NULL,
                             notes   = NULL,
                             label   = NULL,
                             width   = NULL,
                             height  = NULL,
                             path    = NULL) {

    checkmate::assert_string(caption, null.ok = TRUE)
    checkmate::assert_string(notes,   null.ok = TRUE)
    checkmate::assert_string(label,   null.ok = TRUE)
    checkmate::assert_number(width,  lower = 0, upper = 1, null.ok = TRUE)
    checkmate::assert_number(height, lower = 0, upper = 1, null.ok = TRUE)
    checkmate::assert_string(path,    null.ok = TRUE)

    inner_tex <- as.character(tbl)
    inner_tex <- .cm_strip_outer_table_env(inner_tex)

    inner_block <- if (!is.null(path)) {
        dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
        writeLines(inner_tex, path)
        sprintf("\\input{%s}", path)
    } else {
        inner_tex
    }

    sized_block <- if (!is.null(width)) {
        width_spec  <- sprintf("%g\\textwidth", width)
        height_spec <- if (!is.null(height)) sprintf("%g\\textheight", height) else "!"
        sprintf("\\resizebox{%s}{%s}{%s}", width_spec, height_spec, inner_block)
    } else {
        inner_block
    }

    # Wrap the sized block (tabular or resizebox-tabular) in \begin{center}
    # so it is centered, but the centering scope ENDS before the notes line.
    # This way notes flow at the natural left margin without needing
    # \raggedright or a minipage.
    centered_block <- sprintf("\\begin{center}\n%s\n\\end{center}", sized_block)

    cap_line <- if (!is.null(caption)) {
        if (!is.null(label)) {
            sprintf("\\caption{%s}\\label{%s}", caption, label)
        } else {
            sprintf("\\caption{%s}", caption)
        }
    } else ""

    # Force a paragraph break between the table and the notes block.
    notes_line <- if (!is.null(notes)) {
        sprintf(
            "\\par\\vspace{0.5em}\\noindent{\\footnotesize \\textit{Notes}: %s}",
            notes
        )
    } else ""

    tex <- paste(
        "\\begin{table}[H]",
        cap_line,
        centered_block,
        notes_line,
        "\\end{table}",
        sep = "\n"
    )

    structure(
        tex,
        class  = c("cm_latex_block", "knitr_kable", "character"),
        format = "latex"
    )
}

#' Strip an outer table environment if present
#'
#' Removes \begin{table}[...] ... \end{table} wrapping from a LaTeX string
#' if the string is wrapped in one. Also strips any inner \centering and
#' \caption{}/\label{} since wrap_latex_table rebuilds those from
#' user-supplied arguments.
#'
#' @param tex Character. LaTeX source.
#' @return Character. Same string with outer table env removed if present.
#' @keywords internal
.cm_strip_outer_table_env <- function(tex) {
    trimmed <- trimws(tex)

    # (?s) enables dotall — '.' matches across newlines.
    # (.*?) is non-greedy so nested envs don't confuse the match.
    # The placement group [...] is optional and uses (?:) non-capturing.
    pattern <- "(?s)^\\\\begin\\{table\\}(?:\\[[^\\]]*\\])?\\s*(.*?)\\s*\\\\end\\{table\\}$"

    if (grepl(pattern, trimmed, perl = TRUE)) {
        inner <- sub(pattern, "\\1", trimmed, perl = TRUE)

        # Strip a leading \centering (the kable inside already centers itself).
        inner <- sub("^\\s*\\\\centering\\s*\\n?", "", inner, perl = TRUE)

        # Strip any inner \caption{} and optional immediately-following \label{}
        # — wrap_latex_table builds its own from caption/label arguments.
        inner <- gsub(
            "\\\\caption\\{[^}]*\\}(\\\\label\\{[^}]*\\})?\\s*\\n?",
            "",
            inner,
            perl = TRUE
        )

        return(inner)
    }

    tex
}


# R/knit-print.R
# New file. Defines the knit_print methods so chunks render raw LaTeX
# without needing results = "asis".

#' Knitr printing method for causalmetrics LaTeX blocks
#'
#' @param x A `cm_latex_block` object.
#' @param ... Passed to `knitr::asis_output()`.
#'
#' @keywords internal
#' @exportS3Method knitr::knit_print
knit_print.cm_latex_block <- function(x, ...) {
    knitr::asis_output(as.character(x))
}

#' Print method for causalmetrics LaTeX blocks
#'
#' At the R console, `cm_latex_block` objects print as the raw LaTeX code
#' rather than as a character vector wrapped in quotes.
#'
#' @param x A `cm_latex_block` object.
#' @param ... Unused.
#'
#' @keywords internal
#' @export
print.cm_latex_block <- function(x, ...) {
    cat(as.character(x), sep = "\n")
    invisible(x)
}




#' Add a notes line to a kable, in the causalmetrics style
#'
#' Alternative to `kableExtra::footnote()` that uses the simpler
#' `{\\footnotesize \\textit{Notes}: ...}` pattern instead of the
#' `threeparttable` environment. Returns a kable_styling object so it
#' can still chain into other kableExtra functions.
#'
#' Use this when you want notes attached to a kable without going through
#' the full `wrap_latex_table()` wrapping — e.g., for RMarkdown documents
#' that already manage their own table environment.
#'
#' @param tbl A kable or kable_styling object.
#' @param notes Character. The notes content. Will be wrapped in
#'   `{\\footnotesize \\textit{Notes}: ...}`.
#'
#' @return The kable_styling object with a notes line appended.
#'
#' @examples
#' \dontrun{
#' dt <- data.frame(x = 1:3)
#' table_task(dt) |>
#'   kable_notes("Sample of three observations.")
#' }
#'
#' @export
kable_notes <- function(tbl, notes) {
    checkmate::assert_string(notes)
    .cm_check_package("kableExtra")

    current <- as.character(tbl)
    notes_block <- sprintf(
        "\n{\\footnotesize \\textit{Notes}: %s}\n",
        notes
    )

    appended <- paste0(current, notes_block)

    # Re-wrap as kable so downstream functions still recognize it.
    structure(
        appended,
        class  = c("knitr_kable", "character"),
        format = "latex"
    )
}
