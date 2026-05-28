# dev/00_setup_dependencies.R
# Run this once when setting up the development environment.

required <- c(
    "data.table", "fixest", "checkmate", "rlang",
    "generics", "ggplot2"
)

suggested <- c(
    "grf", "rdrobust", "did", "MatchIt",
    "ranger", "xgboost", "glmnet", "estimatr",
    "broom", "modelsummary", "kableExtra", "huxtable",
    "DBI", "knitr", "rmarkdown", "testthat"
)

dev_only <- c(
    "devtools", "usethis", "roxygen2", "pak",
    "renv", "lintr", "styler", "pkgdown", "rcmdcheck"
)

install_if_missing <- function(pkgs) {
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing) > 0L) {
        if (!requireNamespace("pak", quietly = TRUE)) {
            install.packages("pak")
        }
        pak::pkg_install(missing)
    }
    invisible(missing)
}

install_if_missing(required)
install_if_missing(suggested)
install_if_missing(dev_only)
