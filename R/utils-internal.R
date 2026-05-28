# Coerce to data.table (copy if data.frame)
.cm_as_dt <- function(dt) {
    if (data.table::is.data.table(dt)) return(dt)
    if (is.data.frame(dt)) return(data.table::as.data.table(dt))
    rlang::abort("`dt` must be a data.frame or data.table.")
}

# Check required columns exist
.cm_assert_cols <- function(dt, cols) {
    checkmate::assert_data_table(dt)
    missing <- setdiff(cols, names(dt))
    if (length(missing) > 0L) {
        rlang::abort(paste0("Missing columns: ", paste(missing, collapse = ", ")))
    }
}

# Check suggested package is installed
.cm_check_package <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        rlang::abort(
            paste0("Package '", pkg, "' is required. ",
                   "Install with: install.packages('", pkg, "')")
        )
    }
}
