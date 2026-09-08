# Extracted from test-synth.R:204

# prequel ----------------------------------------------------------------------
ferman_sc <- function(y_before, y_after, demean = FALSE) {
  # Ferman and Pinto's replication code (_aux.R), quadprog on raw outcomes
  y <- y_before[, 1]; X <- y_before[, -1]
  if (demean) X <- cbind(1, X)
  Dmat <- t(X) %*% X; dvec <- t(X) %*% y
  if (demean) {
    Amat <- t(rbind(c(0, rep(1, ncol(X) - 1)), cbind(0, diag(ncol(X) - 1)))); bvec <- c(1, rep(0, ncol(X) - 1))
  } else {
    Amat <- t(rbind(rep(1, ncol(X)), diag(ncol(X)))); bvec <- c(1, rep(0, ncol(X)))
  }
  m <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  if (demean) list(w = m$solution[-1], effects = -m$solution[1] + y_after %*% c(1, -m$solution[-1]))
  else list(w = m$solution, effects = y_after %*% c(1, -m$solution))
}

# test -------------------------------------------------------------------------
skip_if_not_installed("quadprog")
skip_if_not_installed("haven")
path <- file.path(testthat::test_path("..", ".."), "data_raw", "09_Synthetic_Control", "andersson2019carbon", "carbontax_data.dta")
skip_if_not(file.exists(path), "Andersson (2019) data not available")
d <- as.data.frame(haven::read_dta(path))
d$treated <- as.integer(d$country == "Sweden" & d$year >= 1990)
set.seed(1)
fit <- synth_control(d, id = "country", time = "year", y = "CO2_transport_capita", d = "treated",
                       x = c("GDP_per_capita", "vehicles_capita", "gas_cons_capita", "urban_pop"),
                       lags = c(1989, 1980, 1970), pre_window = 1980:1989, v = "mspe",
                       v_control = list(maxit = 2000, starts = 3))
published <- c(Denmark = 0.384, Belgium = 0.195, `New Zealand` = 0.177, Greece = 0.090, `United States` = 0.088, Switzerland = 0.061)
w <- setNames(fit$weights$weight, fit$weights$country)
expect_true(all(abs(w[names(published)] - published) < 0.05))
