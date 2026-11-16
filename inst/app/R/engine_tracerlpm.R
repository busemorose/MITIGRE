tlpm_G <- function(type, tau, p2, lam, x) {
  if (type == "PFM") {
    ginf <- exp(-lam * tau)
    return(list(G = ifelse(x >= tau & x > 0, ginf, 0), Ginf = ginf))
  }
  if (type == "DM") {
    sq <- sqrt(1 + 4 * p2 * tau * lam)
    fac <- exp((1 / (2 * p2)) * (1 - sq))
    return(list(G = fac * pinvgauss(x, tau / sq, tau / (2 * p2)), Ginf = fac))
  }
  sr <- switch(type,
               EMM = c(0, 1 / tau),
               EPM = c(tau * (1 - 1 / p2), p2 / tau),
               PEM = {
                 taq <- tau / (1 + log(p2))
                 c(taq * log(p2), 1 / taq)
               })
  s0 <- sr[1]
  r <- sr[2]
  k <- r + lam
  ginf <- r / k * exp(-lam * s0)
  G <- ifelse(x > s0, ginf * (1 - exp(-k * pmax(x - s0, 0))), 0)
  list(G = G, Ginf = ginf)
}

tlpm_kernel_w <- function(type, tau_y, p2, lam, dt, o, uz, nlag) {
  e <- c(0, o * dt + (seq_len(nlag) - 1) * dt)
  g <- tlpm_G(type, tau_y, p2, lam, e - uz)
  G <- g$G
  G[1] <- 0
  f <- exp(-lam * uz)
  list(w = f * diff(G), total = f * g$Ginf, lo = e[-length(e)], hi = e[-1])
}

tlpm_component <- function(type, tau, p2, prep) {
  n <- length(prep$x)
  k <- tlpm_kernel_w(type, tau * prep$dt, p2, prep$lam_year, prep$dt, prep$o, prep$uz, n)
  sim <- apply_kernel(prep$x, k$w, prep$c0, k$total)
  sim[(prep$n_warm + 1):(prep$n_warm + prep$n_obs)]
}

tlpm_kernel <- function(type, tau, p2, prep, decay = FALSE) {
  lam <- if (decay) prep$lam_year else 0
  span <- support_length(type, tau, p2, 1e-7) + prep$uz / prep$dt + prep$o + 2
  nlag <- as.integer(min(ceiling(span), prep$kmax))
  k <- tlpm_kernel_w(type, tau * prep$dt, p2, lam, prep$dt, prep$o, prep$uz, nlag)
  data.frame(lag = seq_len(nlag) - 1L, lo = k$lo, hi = k$hi, w = k$w)
}
