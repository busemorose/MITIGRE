MC_EPM <- function(T, n, isotope, max_t = 1000) {
  t <- 1:max_t
  C <- ifelse(t >= T * (1 - 1 / n), (n / T) * exp(((-n * (t - 0.5)) / T + n - 1)), 0)
  MC <- C / sum(C)
  if (isotope ==  "3h") MC <- MC * exp(-0.0562619464740215 * t)
  return(MC)
}
