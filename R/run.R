validate_comps <- function(comps) {
  for (i in seq_along(comps)) {
    cp <- comps[[i]]
    lab <- paste0("Composante ", i, " (", cp$type, ") : ")
    if (any(is.na(cp$tau)) || cp$tau[1] <= 0 || cp$tau[2] < cp$tau[1]) {
      return(paste0(lab, "il faut 0 < tau min <= tau max."))
    }
    if (cp$type %in% c("EPM", "PEM") && (any(is.na(cp$p2)) || cp$p2[1] < 1 || cp$p2[2] < cp$p2[1])) {
      return(paste0(lab, "il faut 1 <= n min <= n max (n = ratio + 1)."))
    }
    if (cp$type == "DM" && (any(is.na(cp$p2)) || cp$p2[1] <= 0 || cp$p2[2] < cp$p2[1])) {
      return(paste0(lab, "il faut 0 < DP min <= DP max."))
    }
    if (length(comps) > 1 && (any(is.na(cp$ratio)) || cp$ratio[1] < 0 || cp$ratio[2] < cp$ratio[1])) {
      return(paste0(lab, "il faut 0 <= ratio min <= ratio max."))
    }
  }
  NULL
}

run_comparison <- function(engines, comps, prep, obs, eval_idx, crit, allow_NA, n_run, seed,
                           keep_sims_max = 5e6, progress = NULL, cores = 1L) {
  if (.Platform$OS.type == "windows") cores <- 1L
  cores <- max(1L, as.integer(cores))
  par <- sample_params(comps, n_run, seed)
  types <- vapply(comps, function(cp) cp$type, character(1))
  keep_sims <- as.numeric(n_run) * prep$n_obs <= keep_sims_max
  obj <- matrix(NA_real_, n_run, length(engines), dimnames = list(NULL, engines))
  res <- list(engines = engines, types = types, par = par, n_run = n_run, crit = crit,
              allow_NA = allow_NA, eval_idx = eval_idx, prep = prep, keep_sims = keep_sims,
              obj = obj, sims = list(), best = list(),
              time = stats::setNames(numeric(length(engines)), engines), errors = character(0))
  batch <- if (cores > 1L) 10L * cores else max(1L, n_run %/% 50L)
  at_mc <- if (allow_NA) eval_idx[!is.na(obs[eval_idx])] else eval_idx
  for (e in engines) {
    pointwise <- e %in% POINTWISE_ENGINES
    t0 <- proc.time()[["elapsed"]]
    S <- if (keep_sims) matrix(NA_real_, n_run, prep$n_obs) else NULL
    n_err <- 0L
    sim_one <- function(r) {
      tryCatch(simulate_run(e, types, par$tau[r, ], par$p2[r, ], par$ratio[r, ], prep,
                            at = if (pointwise) at_mc else NULL),
               error = function(err) NULL)
    }
    for (bt in split(seq_len(n_run), ceiling(seq_len(n_run) / batch))) {
      sims <- if (cores > 1L) parallel::mclapply(bt, sim_one, mc.cores = cores) else lapply(bt, sim_one)
      for (j in seq_along(bt)) {
        r <- bt[j]
        sim <- sims[[j]]
        if (!is.numeric(sim) || length(sim) != prep$n_obs) {
          n_err <- n_err + 1L
          sim <- rep(NA_real_, prep$n_obs)
        }
        res$obj[r, e] <- crit_value(crit, sim[eval_idx], obs[eval_idx], allow_NA)
        if (keep_sims) S[r, ] <- sim
      }
      if (!is.null(progress)) progress(e, max(bt), length(bt))
    }
    res$time[[e]] <- proc.time()[["elapsed"]] - t0
    if (n_err > 0L) res$errors <- c(res$errors, sprintf("%s : %d run(s) en erreur", ENGINE_LABELS[[e]], n_err))
    res$sims[[e]] <- S
    b <- best_index(res$obj[, e], crit)
    best_sim <- if (is.na(b)) {
      rep(NA_real_, prep$n_obs)
    } else if (keep_sims && !pointwise) {
      S[b, ]
    } else {
      simulate_run(e, types, par$tau[b, ], par$p2[b, ], par$ratio[b, ], prep)
    }
    res$best[[e]] <- list(idx = b, sim = best_sim)
  }
  res
}

run_sim <- function(res, engine, i) {
  if (!is.null(res$sims[[engine]]) && !(engine %in% POINTWISE_ENGINES)) return(res$sims[[engine]][i, ])
  simulate_run(engine, res$types, res$par$tau[i, ], res$par$p2[i, ], res$par$ratio[i, ], res$prep)
}
