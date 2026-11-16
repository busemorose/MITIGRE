mitigre_MC_EMM <- function(T, n, isotope, max_t = 1000) {
  t <- 1:max_t
  C <- (1 / T) * exp(-(t - 0.5) / T)
  MC <- C / sum(C)
  if (isotope == "3h") MC <- MC * exp(-0.0562619464740215 * t)
  return(MC)
}

mitigre_MC_EPM <- function(T, n, isotope, max_t = 1000) {
  t <- 1:max_t
  C <- ifelse(t >= T * (1 - 1 / n), (n / T) * exp(((-n * (t - 0.5)) / T + n - 1)), 0)
  MC <- C / sum(C)
  if (isotope == "3h") MC <- MC * exp(-0.0562619464740215 * t)
  return(MC)
}

mitigre_MC_PEM <- function(T, n, isotope, max_t = 1000) {
  t <- 1:max_t
  C <- ifelse((t - 0.5) >= (T * log(n)), (n / T) * exp(-(t - 0.5) / T), 0)
  MC <- C / sum(C)
  if (isotope == "3h") MC <- MC * exp(-0.0562619464740215 * t)
  return(MC)
}

mitigre_MC_DM <- function(T, DP, isotope, max_t = 1000) {
  t <- 1:max_t
  a <- (t - 0.5) / T
  C <- (1 / ((t - 0.5) * sqrt(4 * pi * DP * a))) * exp(-((1 - a)^2) / (4 * DP * a))
  MC <- C / sum(C)
  if (isotope == "3h") MC <- MC * exp(-0.0562619464740215 * t)
  return(MC)
}

shift_lag <- function(x, k) {
  n <- length(x)
  if (k <= 0) return(x)
  if (k >= n) return(numeric(n))
  c(numeric(k), x[seq_len(n - k)])
}

mitigre_MC_PFM <- function(T, Cin, isotope) {
  T_sup <- ceiling(T)
  sim <- shift_lag(Cin, T_sup)
  if (isotope == "3h") sim <- sim * exp(-0.0562619464740215 * T)
  return(sim)
}

mitigre_native_T <- function(type, tau, p2) {
  if (type == "PEM") tau / (1 + log(p2)) else tau
}

mitigre_MC <- function(type, T, p2, isotope, max_t) {
  switch(type,
         EMM = mitigre_MC_EMM(T, p2, isotope, max_t),
         EPM = mitigre_MC_EPM(T, p2, isotope, max_t),
         PEM = mitigre_MC_PEM(T, p2, isotope, max_t),
         DM = mitigre_MC_DM(T, p2, isotope, max_t))
}

mitigre_component <- function(type, tau, p2, prep) {
  T <- mitigre_native_T(type, tau, p2)
  max_t <- prep$max_t
  if (is.na(max_t)) max_t <- 1
  Cin_warmup <- prep$x
  if (type != "PFM") {
    MC <- mitigre_MC(type, T, p2, prep$isotope, max_t)
    sim <- rev(stats::convolve(rev(MC), Cin_warmup, conj = TRUE, type = "open"))
  } else {
    sim <- mitigre_MC_PFM(T, Cin_warmup, prep$isotope)
  }
  sim[(prep$n_warm + 1):(prep$n_warm + prep$n_obs)]
}

mitigre_kernel <- function(type, tau, p2, prep, decay = FALSE) {
  T <- mitigre_native_T(type, tau, p2)
  iso <- if (decay) prep$isotope else "18o"
  max_t <- prep$max_t
  if (is.na(max_t)) max_t <- 1
  if (type == "PFM") {
    lag <- ceiling(T)
    w <- c(numeric(lag), if (iso == "3h") exp(-0.0562619464740215 * T) else 1)
  } else {
    w <- mitigre_MC(type, T, p2, iso, max_t)
  }
  lag <- seq_along(w) - 1L
  data.frame(lag = lag, lo = lag * prep$dt, hi = (lag + 1) * prep$dt, w = w)
}
