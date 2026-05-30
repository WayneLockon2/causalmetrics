# R/zzz.R

.onLoad <- function(libname, pkgname) {
    # Nothing needed yet; placeholder for future setup
}


# R/zzz.R

#' @keywords internal
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

    # Check which are already attached
    attached <- attach_pkgs %in% .packages()
    to_attach <- attach_pkgs[!attached]

    # Suppress the startup messages from attached packages so the output
    # is clean. The user sees only the causalmetrics banner and a one-line
    # summary of what was attached.
    attach_quietly <- function(pkg) {
        suppressPackageStartupMessages(
            requireNamespace(pkg, quietly = TRUE) &&
                attachNamespace(pkg)
        )
    }

    # We rely on Depends in DESCRIPTION to actually attach these. The block
    # below is defensive — if for some reason a package was unloaded after
    # the Depends mechanism fired, attempt to reattach it.
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

#' @keywords internal
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

#' Detect function name conflicts between attached packages
#' @keywords internal
.cm_detect_conflicts <- function(pkgs) {
    # Returns a character vector of conflict descriptions like
    # "filter(): dplyr masks stats" — but we only check among the
    # packages we attach.

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
