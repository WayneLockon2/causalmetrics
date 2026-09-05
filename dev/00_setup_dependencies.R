# dev/00_setup_dependencies.R
# Run this once when setting up the development environment.

# Runtime dependencies (DESCRIPTION Depends and Imports).
required <- c(
    "tidyverse", "data.table", "fixest", "estimatr", "broom", "modelsummary",
    "kableExtra", "checkmate", "ggplot2", "knitr", "rlang"
)

# Used by est_aipw(), the vignettes, and the tests (DESCRIPTION Suggests).
suggested <- c(
    "bookdown", "car", "dplyr", "ggthemes", "grf", "MatchIt", "mlr3",
    "mlr3learners", "patchwork", "randomizr", "ranger", "ri2", "rmarkdown",
    "sandwich", "testthat", "tibble"
)

# Planned but not yet used: extra mlr3 learners and table/database backends.
planned <- c("glmnet", "xgboost", "huxtable", "DBI")

# Used only by the paper replications in inst/replications/.
replication_only <- c(
    "haven", "scales", "CBPS", "hbal", "Matching", "DoubleML", "qte",
    "sensemakr", "mfx", "DescTools", "ggpubr"
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
install_if_missing(replication_only)
install_if_missing(planned)
install_if_missing(dev_only)
