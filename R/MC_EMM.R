MC_EMM <- function(T, n, isotope, max_t = 1000) {
  # "n" here is useless, it is a lazy workaround for the do.call call in model() function
  t <- 1:max_t
  C <- (1 / T) * exp(-(t - 0.5) / T)
  MC <- C / sum(C)
  if (isotope ==  "3h") MC <- MC * exp(-0.0562619464740215 * t)
  return(MC)
}
