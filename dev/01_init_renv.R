# dev/01_init_renv.R
# Run this once to initialize reproducible environment.

if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv")
}
renv::init()
renv::snapshot()
