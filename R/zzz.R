# R/zzz.R
#
# Package hooks. `.onLoad()` runs whenever the namespace is loaded (including
# `causalmetrics::fun` calls); `.onAttach()` runs only on `library()`.

.onLoad <- function(libname, pkgname) {
    # Use a reasonable default for data.table threading.
    # The user can override with setDTthreads() after loading.
    if (requireNamespace("data.table", quietly = TRUE)) {
        n <- data.table::getDTthreads()
        if (n == 0L) {
            data.table::setDTthreads(percent = 50)
        }
    }
}

.onAttach <- function(libname, pkgname) {
    # Packages to attach when causalmetrics is loaded.
    # These MUST also be listed in DESCRIPTION Depends.
    attach_pkgs <- c(
        "data.table",
        "fixest",
        "estimatr",
        "broom",
        "modelsummary",
        "kableExtra",
        "tidyverse"
    )

    # We rely on Depends in DESCRIPTION to actually attach these. The block
    # below is defensive: if a package was detached after the Depends
    # mechanism fired, attempt to reattach it quietly.
    attached <- attach_pkgs %in% .packages()
    to_attach <- attach_pkgs[!attached]
    for (pkg in to_attach) {
        if (requireNamespace(pkg, quietly = TRUE)) {
            tryCatch(
                suppressPackageStartupMessages(attachNamespace(pkg)),
                error = function(e) invisible(NULL)
            )
        }
    }

    # Startup banner
    ver <- utils::packageVersion("causalmetrics")

    msg <- c(
        sprintf("causalmetrics %s", ver),
        sprintf("Attached: %s", paste(attach_pkgs, collapse = ", "))
    )

    # Surface any conflicts that the user should know about.
    conflicts <- .cm_detect_conflicts(attach_pkgs)
    if (length(conflicts) > 0L) {
        msg <- c(
            msg,
            "",
            "Function name conflicts (later attachments mask earlier ones):",
            paste0("  ", conflicts)
        )
    }

    packageStartupMessage(paste(msg, collapse = "\n"))
}

#' Detect function name conflicts between attached packages
#'
#' @param pkgs Character vector of package names. Only packages currently on
#'   the search path are compared.
#' @return Character vector of conflict descriptions such as
#'   `"filter(): dplyr masks stats"`, or `character(0)` when fewer than two of
#'   the packages are attached or no names overlap.
#' @keywords internal
.cm_detect_conflicts <- function(pkgs) {
    envs <- lapply(pkgs, function(p) {
        nm <- paste0("package:", p)
        if (nm %in% search()) as.environment(nm) else NULL
    })
    names(envs) <- pkgs
    envs <- envs[!vapply(envs, is.null, logical(1))]

    if (length(envs) < 2L) return(character(0))

    # Collect exported function names per package
    exports <- lapply(envs, function(e) ls(e))

    # Find names appearing in more than one package
    all_names <- unlist(exports)
    dup <- unique(all_names[duplicated(all_names)])

    if (length(dup) == 0L) return(character(0))

    # For each duplicated name, report which packages export it
    out <- vapply(dup, function(nm) {
        in_pkgs <- names(exports)[vapply(exports, function(x) nm %in% x, logical(1))]
        sprintf("%s(): %s", nm, paste(rev(in_pkgs), collapse = " masks "))
    }, character(1))

    # Hide some uninteresting conflicts with base R that always appear
    # (these are normal for data.table and tidyverse-style packages)
    boring <- c("first", "last", "between", "transpose")
    out <- out[!names(out) %in% boring]

    unname(out)
}
