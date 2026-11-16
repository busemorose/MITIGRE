tlpm_back_TI <- function(type, tau_y, p2, dt) {
  if (type == "DM") {
    TI <- if (p2 < 1) 1 / 12 else 1 / 24
    if (tau_y >= 100) TI <- dt
  } else {
    TI <- if (tau_y <= 100) 1 / 12 else dt
  }
  m <- max(1, round(dt / TI))
  dt / m
}

tlpm_back_bins <- function(type, tau_y, p2, lam, dt, o, uz, nbin) {
  TI <- tlpm_back_TI(type, tau_y, p2, dt)
  r <- o * dt - floor(o * dt / TI + 1e-9) * TI
  if (r < 1e-12) r <- 0
  e <- c(0, r + (seq_len(nbin) - 1) * TI)
  G <- tlpm_G(type, tau_y, p2, lam, e)$G
  hi <- e[-1]
  list(w = exp(-lam * uz) * diff(G),
       off = pmin(floor((o * dt - (hi + uz)) / dt + 1e-9), 0),
       hi = hi)
}

tlpm_back_component <- function(type, tau, p2, prep, at = NULL) {
  if (is.null(at)) at <- seq_len(prep$n_obs)
  out <- rep(NA_real_, prep$n_obs)
  if (!length(at)) return(out)
  x <- c(prep$c0, prep$x)
  tau_y <- tau * prep$dt
  lam <- prep$lam_year
  pos <- prep$n_warm + at
  if (type == "PFM") {
    off <- min(floor((prep$o * prep$dt - tau_y - prep$uz) / prep$dt + 1e-9), 0)
    kk <- pmax(pos + off + 1, 1)
    out[at] <- x[kk] * exp(-lam * (tau_y + prep$uz))
    return(out)
  }
  fac <- if (type == "DM") 10 else 1
  dt <- prep$dt
  o <- prep$o
  uz <- prep$uz
  TI <- tlpm_back_TI(type, tau_y, p2, dt)
  r <- o * dt - floor(o * dt / TI + 1e-9) * TI
  if (r < 1e-12) r <- 0
  f_uz <- exp(-lam * uz)
  compat <- isTRUE(prep$tlpm_compat) && type == "DM"
  delta <- NULL
  if (compat) {
    TI <- if (p2 < 1) 1 / 12 else 1 / 24
    if (tau_y >= 100) TI <- dt
    S_date <- prep$dates[pos] + o * dt
    End <- trunc(S_date) + floor((S_date - trunc(S_date)) / TI + 1e-9) * TI
    on_grid <- abs(End - S_date) < 1e-9
    End[on_grid] <- End[on_grid] - 1 / 12
    delta <- End - prep$dates[pos]
    r <- 0
  }
  active <- seq_along(at)
  cs <- numeric(length(at))
  G_prev <- if (compat) tlpm_G(type, tau_y, p2, lam, TI)$G else 0
  b0 <- if (compat) 3 else 1
  B <- 256L
  repeat {
    hi <- r + (b0:(b0 + B - 1) - 1) * TI
    G <- tlpm_G(type, tau_y, p2, lam, hi)$G
    w <- f_uz * diff(c(G_prev, G))
    G_prev <- G[B]
    if (compat) {
      idx <- floor(outer(-hi, delta[active], "+") / dt + 1e-9) + rep(pos[active], each = B) + 1
    } else {
      off <- pmin(floor((o * dt - (hi + uz)) / dt + 1e-9), 0)
      idx <- outer(off, pos[active], "+") + 1
    }
    idx[idx < 1] <- 1
    C <- w * matrix(x[idx], nrow = B)
    S <- matrix(cumsum(C), nrow = B)
    if (ncol(S) > 1) S <- S - rep(c(0, S[B, -ncol(S)]), each = B)
    S <- S + rep(cs[active], each = B)
    past <- if (compat) hi > tau_y else hi >= tau_y
    cond <- (abs(S) > 0) & (fac * abs(C) < 1e-6 * abs(S)) & past
    h <- which(cond)
    if (length(h)) {
      hc <- (h - 1L) %/% B + 1L
      first <- !duplicated(hc)
      out[at[active[hc[first]]]] <- S[h[first]]
    }
    cs[active] <- S[B, ]
    if (length(h)) active <- active[-unique(hc)]
    if (!length(active)) break
    b0 <- b0 + B
    if (b0 > 2e6) {
      out[at[active]] <- cs[active]
      break
    }
    B <- min(B * 2L, 8192L)
  }
  out
}

tlpm_back_kernel <- function(type, tau, p2, prep, decay = FALSE) {
  tlpm_kernel(type, tau, p2, prep, decay)
}
