claude_g <- function(type, tau, p2, a) {
  switch(type,
         EMM = (1 / tau) * exp(-a / tau),
         EPM = ifelse(a >= tau * (1 - 1 / p2), (p2 / tau) * exp(-p2 * a / tau + p2 - 1), 0),
         PEM = {
           taq <- tau / (1 + log(p2))
           ifelse(a >= taq * log(p2), (p2 / taq) * exp(-a / taq), 0)
         },
         DM = {
           x <- a / tau
           (1 / (a * sqrt(4 * pi * p2 * x))) * exp(-((1 - x)^2) / (4 * p2 * x))
         })
}

claude_kernel_w <- function(type, tau, p2, lam_step, o, tol, kmax, decay = TRUE) {
  if (!decay) lam_step <- 0
  if (type == "PFM") {
    lag <- max(0, ceiling(tau - o))
    w <- c(numeric(lag), exp(-lam_step * tau))
    lags <- seq_along(w) - 1L
    return(list(w = w, total = sum(w), lag = lags,
                lo = ifelse(lags == 0, 0, o + lags - 1), hi = o + lags, capped = FALSE))
  }
  need <- ceiling(support_length(type, tau, p2, tol) + o) + 1
  K <- as.integer(min(need, kmax))
  lags <- 0:(K - 1L)
  lo <- ifelse(lags == 0, 0, o + lags - 1)
  hi <- o + lags
  mid <- (lo + hi) / 2
  wnd <- claude_g(type, tau, p2, mid) * (hi - lo)
  wnd[!is.finite(wnd)] <- 0
  wnd <- wnd / sum(wnd)
  w <- wnd * exp(-lam_step * mid)
  list(w = w, total = sum(w), lag = lags, lo = lo, hi = hi, capped = need > kmax)
}

claude_component <- function(type, tau, p2, prep) {
  if (identical(prep$quad, "exact")) return(tlpm_component(type, tau, p2, prep))
  k <- claude_kernel_w(type, tau, p2, prep$lam_year * prep$dt, prep$o, prep$tol, prep$kmax)
  sim <- apply_kernel(prep$x, k$w, prep$c0, k$total)
  sim[(prep$n_warm + 1):(prep$n_warm + prep$n_obs)]
}

claude_kernel <- function(type, tau, p2, prep, decay = FALSE) {
  if (identical(prep$quad, "exact")) return(tlpm_kernel(type, tau, p2, prep, decay))
  k <- claude_kernel_w(type, tau, p2, prep$lam_year * prep$dt, prep$o, prep$tol, prep$kmax, decay)
  data.frame(lag = k$lag, lo = k$lo * prep$dt, hi = k$hi * prep$dt, w = k$w)
}
