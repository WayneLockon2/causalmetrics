# R/zzz.R

.onLoad <- function(libname, pkgname) {
    # Nothing needed yet; placeholder for future setup
}

.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        "causalmetrics ", utils::packageVersion("causalmetrics"),
        " | Backend: data.table + fixest"
    )
}
