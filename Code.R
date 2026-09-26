# ---------- 0) Utilities ----------

expit <- function(z) 1/(1+exp(-z))

mu_fun <- function(eta, family) {
  if (family == "logistic") return(expit(eta))
  if (family == "poisson")  return(exp(eta))
  stop("family must be 'logistic' or 'poisson'")
}

dmu_fun <- function(eta, family) {
  if (family == "logistic") {
    m <- expit(eta)
    return(m * (1 - m))
  }
  if (family == "poisson") return(exp(eta))
  stop("family must be 'logistic' or 'poisson'")
}

summarize_mc <- function(est_mat, true_psi) {
  # est_mat: reps x p matrix
  # true_psi: length p
  stopifnot(ncol(est_mat) == length(true_psi))
  
  bias <- colMeans(est_mat - matrix(true_psi, nrow(est_mat), length(true_psi), byrow = TRUE), na.rm = TRUE)
  sd   <- apply(est_mat, 2, sd, na.rm = TRUE)
  mse  <- colMeans((est_mat - matrix(true_psi, nrow(est_mat), length(true_psi), byrow = TRUE))^2, na.rm = TRUE)
  
  data.frame(param = names(true_psi), truth = true_psi, bias = bias, sd = sd, mse = mse)
}

summarize_mc_with_ci <- function(est_mat, se_mat, true_psi) {
  stopifnot(ncol(est_mat) == length(true_psi))
  stopifnot(all(dim(est_mat) == dim(se_mat)))
  
  truth_mat <- matrix(true_psi, nrow(est_mat), length(true_psi), byrow = TRUE)
  cover_mat <- (est_mat - 1.96 * se_mat <= truth_mat) & (truth_mat <= est_mat + 1.96 * se_mat)
  
  data.frame(
    param = names(true_psi),
    truth = true_psi,
    bias = colMeans(est_mat - truth_mat, na.rm = TRUE),
    sd = apply(est_mat, 2, sd, na.rm = TRUE),
    mean_se = colMeans(se_mat, na.rm = TRUE),
    coverage = colMeans(cover_mat, na.rm = TRUE),
    mse = colMeans((est_mat - truth_mat)^2, na.rm = TRUE)
  )
}

summarize_scalar_mc_with_ci <- function(est, se, truth) {
  cover <- (est - 1.96 * se <= truth) & (truth <= est + 1.96 * se)
  c(
    truth = mean(truth, na.rm = TRUE),
    bias = mean(est - truth, na.rm = TRUE),
    sd = sd(est, na.rm = TRUE),
    mean_se = mean(se, na.rm = TRUE),
    coverage = mean(cover, na.rm = TRUE),
    mse = mean((est - truth)^2, na.rm = TRUE)
  )
}

as_family_object <- function(family) {
  if (family == "logistic") return(binomial())
  if (family == "poisson") return(poisson())
  stop("family must be 'logistic' or 'poisson'")
}

binary_Pi <- function(p01, p10) {
  matrix(c(1 - p01, p01,
           p10, 1 - p10), nrow = 2, ncol = 2)
}

validate_mismeasured <- function() {
  if (!requireNamespace("mismeasured", quietly = TRUE)) {
    stop("Package 'mismeasured' is required. Install it with remotes::install_github('OJALAB/mismeasured').")
  }
}

safe_vcov <- function(object, method = NULL) {
  out <- tryCatch(vcov(object, method = method), error = function(e) NULL)
  if (is.null(out)) return(NULL)
  out
}

reorder_simex_binary <- function(beta, V) {
  beta_out <- c(gamma = unname(beta["1"]),
                alpha0 = unname(beta["(Intercept)"]),
                alpha1 = unname(beta["x1"]))
  V_out <- matrix(NA_real_, 3, 3, dimnames = list(names(beta_out), names(beta_out)))
  if (!is.null(V)) {
    old <- c("1", "(Intercept)", "x1")
    V_out[,] <- V[old, old, drop = FALSE]
  }
  list(coef = beta_out, vcov = V_out)
}

reorder_simex_K <- function(beta, V, K) {
  old <- c(as.character(1:(K - 1)), "(Intercept)", "x1")
  new <- c(paste0("gamma", 1:(K - 1)), "alpha0", "alpha1")
  beta_out <- setNames(as.numeric(beta[old]), new)
  V_out <- matrix(NA_real_, length(new), length(new), dimnames = list(new, new))
  if (!is.null(V)) V_out[,] <- V[old, old, drop = FALSE]
  list(coef = beta_out, vcov = V_out)
}

first_existing_name <- function(candidates, available) {
  hit <- candidates[candidates %in% available]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

standardize_mcglm_binary <- function(beta, V = NULL) {
  old <- c(
    gamma = first_existing_name(c("gamma", "1"), names(beta)),
    alpha0 = first_existing_name(c("alpha0", "(Intercept)", "Intercept"), names(beta)),
    alpha1 = first_existing_name(c("alpha1", "x1"), names(beta))
  )
  new <- names(old)
  
  beta_out <- setNames(rep(NA_real_, length(new)), new)
  ok <- !is.na(old)
  beta_out[ok] <- as.numeric(beta[old[ok]])
  
  V_out <- matrix(NA_real_, length(new), length(new), dimnames = list(new, new))
  if (!is.null(V)) {
    ok_v <- ok & old %in% rownames(V) & old %in% colnames(V)
    V_out[ok_v, ok_v] <- V[old[ok_v], old[ok_v], drop = FALSE]
  }
  
  list(coef = beta_out, vcov = V_out)
}

standardize_mcglm_K <- function(beta, V = NULL, K) {
  gamma_old <- vapply(1:(K - 1), function(k) {
    first_existing_name(c(paste0("gamma", k), paste0("gamma_", k), as.character(k)), names(beta))
  }, character(1))
  old <- c(
    setNames(gamma_old, paste0("gamma", 1:(K - 1))),
    alpha0 = first_existing_name(c("alpha0", "(Intercept)", "Intercept"), names(beta)),
    alpha1 = first_existing_name(c("alpha1", "x1"), names(beta))
  )
  new <- names(old)
  
  beta_out <- setNames(rep(NA_real_, length(new)), new)
  ok <- !is.na(old)
  beta_out[ok] <- as.numeric(beta[old[ok]])
  
  V_out <- matrix(NA_real_, length(new), length(new), dimnames = list(new, new))
  if (!is.null(V)) {
    ok_v <- ok & old %in% rownames(V) & old %in% colnames(V)
    V_out[ok_v, ok_v] <- V[old[ok_v], old[ok_v], drop = FALSE]
  }
  
  list(coef = beta_out, vcov = V_out)
}

fit_methods_binary_pkg <- function(Y, W, X, family, p01, p10, pi,
                                   include_simex = TRUE, simex_B = 25,
                                   simex_method = "improved", simex_seed = 1, validation = NULL) {
  validate_mismeasured()
  mc <- mismeasured::mc
  
  Pi <- binary_Pi(p01, p10)
  dat <- data.frame(Y = Y, W = factor(W, levels = 0:1), x1 = X[, 2])
  methods <- c("naive", "bca", "bcm", "cs", "cs_akn")
  pnames <- c("gamma", "alpha0", "alpha1")
  
  coefs <- matrix(NA_real_, nrow = length(methods), ncol = length(pnames),
                  dimnames = list(methods, pnames))
  vcovs <- setNames(vector("list", length(methods)), methods)
  if(is.null(validation)){
    fit <- tryCatch(
      mismeasured::mcglm(Y ~ mc(W, Pi) + x1, data = dat, family = as_family_object(family),
                         method = methods, pi_z = pi, vcov_corrected = TRUE),
      error = function(e) NULL
    )
  } else{
    fit <- mismeasured::mcglm(Y, z_hat = W, x = X, family = as_family_object(family), 
                       method = methods, Pi = Pi, pi_z = pi, vcov_corrected = TRUE, validation = validation)
  }
  if (!is.null(fit)) {
    for (m in methods) {
      std <- standardize_mcglm_binary(fit$coefficients[[m]], safe_vcov(fit, method = m))
      coefs[m, ] <- std$coef[pnames]
      vcovs[[m]] <- std$vcov[pnames, pnames, drop = FALSE]
    }
  }
  
  if (include_simex) {
    methods <- c(methods, "simex")
    coefs <- rbind(coefs, simex = setNames(rep(NA_real_, length(pnames)), pnames))
    vcovs[["simex"]] <- matrix(NA_real_, length(pnames), length(pnames), dimnames = list(pnames, pnames))
    
    fit_simex <- tryCatch(
      mismeasured::simex(Y ~ mc(W, Pi) + x1, family = as_family_object(family),
                         data = dat, method = simex_method, B = simex_B,
                         jackknife = TRUE, seed = simex_seed),
      error = function(e) NULL
    )
    if (!is.null(fit_simex)) {
      sim <- reorder_simex_binary(coef(fit_simex), safe_vcov(fit_simex))
      coefs["simex", ] <- sim$coef
      vcovs[["simex"]] <- sim$vcov
    }
  }
  
  list(coef = coefs, vcov = vcovs)
}

fit_methods_K_pkg <- function(Y, W, X, K, family, Pi, pi_vec,
                              include_simex = TRUE, simex_B = 25,
                              simex_method = "improved", simex_seed = 1) {
  validate_mismeasured()
  mc <- mismeasured::mc
  
  dat <- data.frame(Y = Y, W = factor(W, levels = 0:(K - 1)), x1 = X[, 2])
  methods <- c("naive", "bca", "bcm", "cs", "cs_akn")
  pnames <- c(paste0("gamma", 1:(K - 1)), "alpha0", "alpha1")
  
  coefs <- matrix(NA_real_, nrow = length(methods), ncol = length(pnames),
                  dimnames = list(methods, pnames))
  vcovs <- setNames(vector("list", length(methods)), methods)
  
  fit <- tryCatch(
    mismeasured::mcglm(Y ~ mc(W, Pi) + x1, data = dat, family = as_family_object(family),
                       method = methods, pi_z = pi_vec, vcov_corrected = TRUE),
    error = function(e) NULL
  )
  
  if (!is.null(fit)) {
    for (m in methods) {
      std <- standardize_mcglm_K(fit$coefficients[[m]], safe_vcov(fit, method = m), K)
      coefs[m, ] <- std$coef[pnames]
      vcovs[[m]] <- std$vcov[pnames, pnames, drop = FALSE]
    }
  }
  
  if (include_simex) {
    methods <- c(methods, "simex")
    coefs <- rbind(coefs, simex = setNames(rep(NA_real_, length(pnames)), pnames))
    vcovs[["simex"]] <- matrix(NA_real_, length(pnames), length(pnames), dimnames = list(pnames, pnames))
    
    fit_simex <- tryCatch(
      mismeasured::simex(Y ~ mc(W, Pi) + x1, family = as_family_object(family),
                         data = dat, method = simex_method, B = simex_B,
                         jackknife = TRUE, seed = simex_seed),
      error = function(e) NULL
    )
    if (!is.null(fit_simex)) {
      sim <- reorder_simex_K(coef(fit_simex), safe_vcov(fit_simex), K)
      coefs["simex", ] <- sim$coef
      vcovs[["simex"]] <- sim$vcov
    }
  }
  
  list(coef = coefs, vcov = vcovs)
}

# ---------- 1) Data generators ----------

simulate_binary <- function(n, family, alpha, gamma, pi, p01, p10) {
  # X: intercept + one covariate
  x1 <- rnorm(n)
  X  <- cbind(1, x1)
  
  # latent Z independent of X
  Z  <- rbinom(n, 1, pi)
  
  # outcome
  eta <- drop(X %*% alpha + gamma * Z)
  mu  <- mu_fun(eta, family)
  Y   <- if (family == "logistic") rbinom(n, 1, mu) else rpois(n, mu)
  
  # misclassify Z -> W
  W <- Z
  idx0 <- which(Z == 0)
  idx1 <- which(Z == 1)
  if (length(idx0) > 0) W[idx0] <- rbinom(length(idx0), 1, p01)         # W=1 with prob p01, else 0
  if (length(idx1) > 0) W[idx1] <- rbinom(length(idx1), 1, 1 - p10)     # W=1 with prob 1-p10, else 0
  
  list(Y = Y, X = X, Z = Z, W = W)
}

estimate_binary_misclass <- function(Z, W) {
  # Estimates p01 = P(W=1|Z=0), p10 = P(W=0|Z=1), and pi = P(Z=1)
  p01_hat <- if (sum(Z == 0) > 0) mean(W[Z == 0] == 1) else NA_real_
  p10_hat <- if (sum(Z == 1) > 0) mean(W[Z == 1] == 0) else NA_real_
  pi_hat  <- mean(Z == 1)
  list(p01 = p01_hat, p10 = p10_hat, pi = pi_hat)
}

simulate_K <- function(n, K, family, alpha, gamma_vec, pi_vec, Pi) {
  # Z takes values 0,...,K-1. Baseline category is 0.
  stopifnot(length(pi_vec) == K)
  stopifnot(length(gamma_vec) == K - 1)
  stopifnot(all(dim(Pi) == c(K, K)))
  
  # X: intercept + one covariate
  x1 <- rnorm(n)
  X  <- cbind(1, x1)
  
  # latent Z independent of X
  Z <- sample(0:(K - 1), size = n, replace = TRUE, prob = pi_vec)
  
  # outcome
  gamma_full <- c(0, gamma_vec)  # gamma_0 = 0
  eta <- drop(X %*% alpha + gamma_full[Z + 1])
  mu  <- mu_fun(eta, family)
  Y   <- if (family == "logistic") rbinom(n, 1, mu) else rpois(n, mu)
  
  # misclassify: W|Z=l ~ Multinomial(Pi[,l])
  W <- vapply(Z, function(z) sample(0:(K - 1), size = 1, prob = Pi[, z + 1]), numeric(1))
  
  list(Y = Y, X = X, Z = Z, W = W)
}

estimate_K_misclass <- function(Z, W, K) {
  # Estimates Pi_{j,l} = P(W=j|Z=l) and class prevalences pi_l = P(Z=l)
  Pi_hat <- matrix(0, nrow = K, ncol = K)
  for (l in 0:(K - 1)) {
    idx <- which(Z == l)
    if (length(idx) == 0) {
      Pi_hat[, l + 1] <- NA_real_
    } else {
      tab <- table(factor(W[idx], levels = 0:(K - 1)))
      Pi_hat[, l + 1] <- as.numeric(tab) / length(idx)
    }
  }
  pi_hat <- as.numeric(table(factor(Z, levels = 0:(K - 1)))) / length(Z)
  list(Pi = Pi_hat, pi = pi_hat)
}

# ---------- 2) Monte Carlo drivers ----------

run_mc_binary <- function(reps = 100, n = 10000, family = "logistic",
                          alpha = c(-0.5, 0.7), gamma = 0.8,
                          pi = 0.4, p01 = 0.10, p10 = 0.15,
                          known_probs = TRUE, val_frac = 0.25,
                          drifting = FALSE, drift_scale = c(2, 2),
                          include_simex = TRUE, simex_B = 25,
                          simex_method = "improved",
                          seed = 1) {
  set.seed(seed)
  
  # If drifting = TRUE, use p01 = c01/sqrt(n) and p10 = c10/sqrt(n)
  if (drifting) {
    p01 <- drift_scale[1] / sqrt(n)
    p10 <- drift_scale[2] / sqrt(n)
  }
  
  true_psi <- c(gamma = gamma, alpha0 = alpha[1], alpha1 = alpha[2])
  
  # store estimates: array reps x methods x p
  methods <- c("naive", "bca", "bcm", "cs", "cs_akn")
  if (include_simex) methods <- c(methods, "simex")
  est_arr <- array(NA_real_, dim = c(reps, length(methods), length(true_psi)),
                   dimnames = list(NULL, methods, names(true_psi)))
  se_arr <- est_arr
  
  for (b in 1:reps) {
    dat <- simulate_binary(n, family, alpha, gamma, pi, p01, p10)
    
    # Choose misclassification params (oracle or plug-in from validation data)
    if (known_probs) {
      p01_use <- p01; p10_use <- p10; pi_use <- pi
      fit <- fit_methods_binary_pkg(dat$Y, dat$W, dat$X, family, p01_use, p10_use, pi_use,
                                    include_simex = include_simex, simex_B = simex_B,
                                    simex_method = simex_method, simex_seed = seed + b)
    } else {
      # Validation subsample where true Z is observed.
      # With val_frac fixed, p-hats are sqrt(n)-consistent: O_p(n^{-1/2}), which is stronger than o_p(n^{-1/4}).
      idx_val <- sample.int(n, size = max(50, floor(val_frac * n)))
      est <- estimate_binary_misclass(dat$Z[idx_val], dat$W[idx_val])
      p01_use <- est$p01; p10_use <- est$p10; pi_use <- est$pi
      
      # protect against degenerate validation splits
      if (!is.finite(p01_use) || !is.finite(p10_use) || !is.finite(pi_use)) {
        next
      }
      fit <- fit_methods_binary_pkg(dat$Y, dat$W, dat$X, family, p01_use, p10_use, pi_use,
                                    include_simex = include_simex, simex_B = simex_B,
                                    simex_method = simex_method, simex_seed = seed + b,
                                    validation = list(z = dat$Z[idx_val], index = idx_val))
    }
    
    est_arr[b, , ] <- fit$coef[methods, names(true_psi), drop = FALSE]
    for (m in methods) {
      V <- fit$vcov[[m]]
      if (!is.null(V) && all(dim(V) == c(length(true_psi), length(true_psi)))) {
        se_arr[b, m, ] <- sqrt(pmax(diag(V[names(true_psi), names(true_psi), drop = FALSE]), 0))
      }
    }
  }
  
  truth_arr <- array(rep(true_psi, each = reps * length(methods)),
                     dim = dim(est_arr), dimnames = dimnames(est_arr))
  cover_arr <- (est_arr - 1.96 * se_arr <= truth_arr) & (truth_arr <= est_arr + 1.96 * se_arr)
  
  # Return raw and summaries
  list(
    settings = list(reps = reps, n = n, family = family, alpha = alpha, gamma = gamma,
                    pi = pi, p01 = p01, p10 = p10, known_probs = known_probs,
                    drifting = drifting, include_simex = include_simex,
                    simex_B = simex_B, simex_method = simex_method),
    true_psi = true_psi,
    est_arr  = est_arr,
    se_arr = se_arr,
    cover_arr = cover_arr
  )
}

run_mc_K <- function(reps = 100, n = 10000, K = 3, family = "logistic",
                     alpha = c(-0.3, 0.5), gamma_vec = c(0.7, -0.4),
                     pi_vec = NULL, Pi = NULL,
                     known_probs = TRUE, val_frac = 0.25,
                     drifting = FALSE, drift_scale = 2,
                     include_simex = TRUE, simex_B = 25,
                     simex_method = "improved",
                     seed = 1) {
  set.seed(seed)
  
  if (is.null(pi_vec)) pi_vec <- rep(1 / K, K)
  
  if (is.null(Pi)) {
    # default: symmetric misclassification
    # P(correct) = 1 - p, P(misclass to any other) = p/(K-1)
    p <- 0.20
    if (drifting) p <- drift_scale / sqrt(n)
    Pi <- matrix(p / (K - 1), nrow = K, ncol = K)
    diag(Pi) <- 1 - p
  }
  
  true_psi <- c(setNames(gamma_vec, paste0("gamma", 1:(K - 1))),
                alpha0 = alpha[1], alpha1 = alpha[2])
  
  methods <- c("naive", "bca", "bcm", "cs", "cs_akn")
  if (include_simex) methods <- c(methods, "simex")
  est_arr <- array(NA_real_, dim = c(reps, length(methods), length(true_psi)),
                   dimnames = list(NULL, methods, names(true_psi)))
  se_arr <- est_arr
  
  for (b in 1:reps) {
    dat <- simulate_K(n, K, family, alpha, gamma_vec, pi_vec, Pi)
    
    if (known_probs) {
      Pi_use <- Pi; pi_use <- pi_vec
    } else {
      idx_val <- sample.int(n, size = max(80, floor(val_frac * n)))
      est <- estimate_K_misclass(dat$Z[idx_val], dat$W[idx_val], K)
      Pi_use <- est$Pi; pi_use <- est$pi
      if (any(!is.finite(Pi_use)) || any(!is.finite(pi_use))) next
    }
    
    fit <- fit_methods_K_pkg(dat$Y, dat$W, dat$X, K, family, Pi_use, pi_use,
                             include_simex = include_simex, simex_B = simex_B,
                             simex_method = simex_method, simex_seed = seed + b)
    est_arr[b, , ] <- fit$coef[methods, names(true_psi), drop = FALSE]
    for (m in methods) {
      V <- fit$vcov[[m]]
      if (!is.null(V) && all(dim(V) == c(length(true_psi), length(true_psi)))) {
        se_arr[b, m, ] <- sqrt(pmax(diag(V[names(true_psi), names(true_psi), drop = FALSE]), 0))
      }
    }
  }
  
  truth_arr <- array(rep(true_psi, each = reps * length(methods)),
                     dim = dim(est_arr), dimnames = dimnames(est_arr))
  cover_arr <- (est_arr - 1.96 * se_arr <= truth_arr) & (truth_arr <= est_arr + 1.96 * se_arr)
  
  list(
    settings = list(reps = reps, n = n, K = K, family = family, alpha = alpha, gamma_vec = gamma_vec,
                    pi_vec = pi_vec, Pi = Pi, known_probs = known_probs, drifting = drifting,
                    include_simex = include_simex, simex_B = simex_B,
                    simex_method = simex_method),
    true_psi = true_psi,
    est_arr  = est_arr,
    se_arr = se_arr,
    cover_arr = cover_arr
  )
}

# ---------- 3) Convenience: print compact summaries ----------

print_method_summaries <- function(mc_obj) {
  est_arr <- mc_obj$est_arr
  true_psi <- mc_obj$true_psi
  
  methods <- dimnames(est_arr)[[2]]
  pnames  <- dimnames(est_arr)[[3]]
  
  cat("\n=== Monte Carlo summary ===\n")
  print(mc_obj$settings)
  cat("\nTrue parameters:\n")
  print(true_psi)
  
  for (m in methods) {
    est_mat <- est_arr[, m, , drop = FALSE]
    est_mat <- matrix(est_mat, ncol = length(pnames))
    colnames(est_mat) <- pnames
    
    cat("\n---", m, "---\n")
    if (!is.null(mc_obj$se_arr)) {
      se_mat <- mc_obj$se_arr[, m, , drop = FALSE]
      se_mat <- matrix(se_mat, ncol = length(pnames))
      colnames(se_mat) <- pnames
      print(summarize_mc_with_ci(est_mat, se_mat, true_psi), row.names = FALSE)
    } else {
      print(summarize_mc(est_mat, true_psi), row.names = FALSE)
    }
  }
  
  invisible(NULL)
}

# ---------- 4) Combining probability and nonprobability samples ----------

draw_nonprob_indices_by_x <- function(x, n_nonprob, nonprob_ratio) {
  n_large <- round(n_nonprob * nonprob_ratio)
  n_small <- n_nonprob - n_large
  
  idx_large <- which(x > 0)
  idx_small <- which(x <= 0)
  
  if (n_large > length(idx_large) || n_small > length(idx_small)) {
    stop("Requested nonprobability sample is larger than the available X strata.")
  }
  
  c(sample(idx_small, n_small), sample(idx_large, n_large))
}

predict_binary_true_Z_mean <- function(Z, X, params, family = "logistic") {
  gamma <- unname(params["gamma"])
  alpha <- unname(params[c("alpha0", "alpha1")])
  
  eta <- gamma * Z + drop(X %*% alpha)
  mu_fun(eta, family)
}

predict_K_true_Z_mean <- function(Z, X, params, K, family = "logistic") {
  r <- ncol(X)
  gamma_full <- c(0, unname(params[1:(K - 1)]))
  alpha <- unname(params[K:(K - 1 + r)])
  
  eta <- drop(X %*% alpha) + gamma_full[Z + 1]
  mu_fun(eta, family)
}

mass_variance_binary <- function(Z, X, params, V, family, n_population) {
  pred <- predict_binary_true_Z_mean(Z, X, params, family)
  n_prob <- length(pred)
  fpc <- 1 - n_prob / n_population
  Vp <- if (n_prob > 1) fpc * stats::var(pred) / n_prob else 0
  
  if (is.null(V) || any(!is.finite(V))) return(Vp)
  xi <- cbind(Z, X)
  colnames(xi) <- c("gamma", "alpha0", "alpha1")
  eta <- drop(xi %*% params[colnames(xi)])
  G <- colMeans(xi * dmu_fun(eta, family))
  Vnp <- as.numeric(t(G) %*% V[colnames(xi), colnames(xi), drop = FALSE] %*% G)
  max(Vp + Vnp, 0)
}

mass_variance_K <- function(Z, X, params, V, K, family, n_population) {
  pred <- predict_K_true_Z_mean(Z, X, params, K, family)
  n_prob <- length(pred)
  fpc <- 1 - n_prob / n_population
  Vp <- if (n_prob > 1) fpc * stats::var(pred) / n_prob else 0
  
  if (is.null(V) || any(!is.finite(V))) return(Vp)
  D <- sapply(1:(K - 1), function(k) as.numeric(Z == k))
  if (K - 1 == 1) D <- matrix(D, ncol = 1)
  xi <- cbind(D, X)
  colnames(xi) <- c(paste0("gamma", 1:(K - 1)), "alpha0", "alpha1")
  eta <- drop(xi %*% params[colnames(xi)])
  G <- colMeans(xi * dmu_fun(eta, family))
  Vnp <- as.numeric(t(G) %*% V[colnames(xi), colnames(xi), drop = FALSE] %*% G)
  max(Vp + Vnp, 0)
}

gold_mean_variance <- function(Y_sample, n_population) {
  n_prob <- length(Y_sample)
  fpc <- 1 - n_prob / n_population
  if (n_prob > 1) fpc * stats::var(Y_sample) / n_prob else 0
}

naive_nonprob_mean_variance <- function(Y_sample) {
  n_sample <- length(Y_sample)
  if (n_sample > 1) stats::var(Y_sample) / n_sample else 0
}

run_comb_binary <- function(reps = 100, n = 10000, n_nonprob = 5000, n_prob = 1000, nonprob_ratio = 0.7,
                            family = 'logistic', alpha = c(-0.5, 0.7), gamma = 0.8,
                            pi = 0.4, p01 = 0.1, p10 = 0.15,
                            known_probs = FALSE, val_frac = 0.1,
                            include_simex = TRUE, simex_B = 25,
                            simex_method = "improved", seed = 1){
  set.seed(seed)
  
  # store estimates: array reps x methods x p
  methods <- c("gold", "naive", "naive MI", "bca", "bcm", "cs", "cs_akn")
  if (include_simex) methods <- c(methods, "simex")
  mean_est <- matrix(NA_real_, nrow = reps, ncol = length(methods),
                     dimnames = list(NULL, methods))
  mean_se <- mean_est
  truth <- rep(NA_real_, reps)
  
  for (b in 1:reps) {
    population <- simulate_binary(n, family, alpha, gamma, pi, p01, p10)
    prob_idx <- sample.int(n, n_prob)
    nonprob_idx <- draw_nonprob_indices_by_x(population$X[, 2], n_nonprob, nonprob_ratio)
    truth[b] <- mean(population$Y)
    
    # Choose correction parameters.
    if (known_probs) {
      p01_use <- p01; p10_use <- p10; pi_use <- pi
      validation <- NULL
    } else {
      val_size <- min(length(nonprob_idx), max(50, floor(val_frac * length(nonprob_idx))))
      idx_val <- sample(n_nonprob, size = val_size)
      est <- estimate_binary_misclass(population$Z[nonprob_idx][idx_val], population$W[nonprob_idx][idx_val])
      p01_use <- est$p01; p10_use <- est$p10; pi_use <- est$pi
      validation <- list(z = population$Z[nonprob_idx][idx_val], index = idx_val)
    }
    
    if (!is.finite(p01_use) || !is.finite(p10_use) || !is.finite(pi_use)) next
    
    fit <- fit_methods_binary_pkg(population$Y[nonprob_idx],
                                  population$W[nonprob_idx],
                                  population$X[nonprob_idx, , drop = FALSE],
                                  family, p01_use, p10_use, pi_use,
                                  include_simex = include_simex, simex_B = simex_B,
                                  simex_method = simex_method, simex_seed = seed + b, validation = validation)
    
    for (i in methods) {
      if (i == "gold") {
        mean_est[b, i] <- mean(population$Y[prob_idx])
        mean_se[b, i] <- sqrt(gold_mean_variance(population$Y[prob_idx], n))
      } else if (i == "naive") {
        mean_est[b, i] <- mean(population$Y[nonprob_idx])
        mean_se[b, i] <- sqrt(naive_nonprob_mean_variance(population$Y[nonprob_idx]))
      } else {
        fit_method <- if (i == "naive MI") "naive" else i
        params <- fit$coef[fit_method, ]
        mean_est[b, i] <- mean(predict_binary_true_Z_mean(population$Z[prob_idx],
                                                          population$X[prob_idx, , drop = FALSE],
                                                          params, family))
        mean_se[b, i] <- sqrt(mass_variance_binary(population$Z[prob_idx],
                                                   population$X[prob_idx, , drop = FALSE],
                                                   params, fit$vcov[[fit_method]], family, n))
      }
    }
  }
  
  cover <- (mean_est - 1.96 * mean_se <= truth) & (truth <= mean_est + 1.96 * mean_se)
  errors <- array(mean_est - truth, dim = c(reps, length(methods), 1),
                  dimnames = list(NULL, methods, "error"))
  
  # Return raw and summaries
  list(
    settings = list(reps = reps, n = n, n_prob = n_prob, n_nonprob = n_nonprob, family = family,
                    alpha = alpha, gamma = gamma, pi = pi, p01 = p01, p10 = p10,
                    nonprob_ratio = nonprob_ratio, known_probs = known_probs, val_frac = val_frac,
                    include_simex = include_simex, simex_B = simex_B,
                    simex_method = simex_method),
    true_mean = truth,
    mean_est = mean_est,
    mean_se = mean_se,
    cover = cover,
    errors  = errors
  )
}

run_comb_K <- function(reps = 100, n = 10000, n_nonprob = 5000, n_prob = 1000, nonprob_ratio = 0.7,
                       K = 3, family = "logistic",
                       alpha = c(-0.3, 0.5), gamma_vec = c(0.7, -0.4),
                       pi_vec = NULL, Pi = NULL,
                       known_probs = FALSE, val_frac = 0.1,
                       include_simex = TRUE, simex_B = 25,
                       simex_method = "improved", seed = 1) {
  set.seed(seed)
  
  if (is.null(pi_vec)) pi_vec <- rep(1 / K, K)
  
  if (is.null(Pi)) {
    p <- 0.20
    Pi <- matrix(p / (K - 1), nrow = K, ncol = K)
    diag(Pi) <- 1 - p
  }
  
  stopifnot(length(gamma_vec) == K - 1)
  stopifnot(length(pi_vec) == K)
  stopifnot(all(dim(Pi) == c(K, K)))
  
  methods <- c("gold", "naive", "naive MI", "bca", "bcm", "cs", "cs_akn")
  if (include_simex) methods <- c(methods, "simex")
  mean_est <- matrix(NA_real_, nrow = reps, ncol = length(methods),
                     dimnames = list(NULL, methods))
  mean_se <- mean_est
  truth <- rep(NA_real_, reps)
  
  for (b in 1:reps) {
    population <- simulate_K(n, K, family, alpha, gamma_vec, pi_vec, Pi)
    prob_idx <- sample.int(n, n_prob)
    nonprob_idx <- draw_nonprob_indices_by_x(population$X[, 2], n_nonprob, nonprob_ratio)
    truth[b] <- mean(population$Y)
    
    if (known_probs) {
      Pi_use <- Pi
      pi_use <- pi_vec
    } else {
      val_size <- min(length(nonprob_idx), max(80, floor(val_frac * length(nonprob_idx))))
      idx_val <- sample(nonprob_idx, size = val_size)
      est <- estimate_K_misclass(population$Z[idx_val], population$W[idx_val], K)
      Pi_use <- est$Pi
      pi_use <- est$pi
    }
    
    if (any(!is.finite(Pi_use)) || any(!is.finite(pi_use))) next
    
    fit <- fit_methods_K_pkg(population$Y[nonprob_idx],
                             population$W[nonprob_idx],
                             population$X[nonprob_idx, , drop = FALSE],
                             K, family, Pi_use, pi_use,
                             include_simex = include_simex, simex_B = simex_B,
                             simex_method = simex_method, simex_seed = seed + b)
    
    for (i in methods) {
      if (i == "gold") {
        mean_est[b, i] <- mean(population$Y[prob_idx])
        mean_se[b, i] <- sqrt(gold_mean_variance(population$Y[prob_idx], n))
      } else if (i == "naive") {
        mean_est[b, i] <- mean(population$Y[nonprob_idx])
        mean_se[b, i] <- sqrt(naive_nonprob_mean_variance(population$Y[nonprob_idx]))
      } else {
        fit_method <- if (i == "naive MI") "naive" else i
        params <- fit$coef[fit_method, ]
        mean_est[b, i] <- mean(predict_K_true_Z_mean(population$Z[prob_idx],
                                                     population$X[prob_idx, , drop = FALSE],
                                                     params, K, family))
        mean_se[b, i] <- sqrt(mass_variance_K(population$Z[prob_idx],
                                              population$X[prob_idx, , drop = FALSE],
                                              params, fit$vcov[[fit_method]], K, family, n))
      }
    }
  }
  
  cover <- (mean_est - 1.96 * mean_se <= truth) & (truth <= mean_est + 1.96 * mean_se)
  errors <- array(mean_est - truth, dim = c(reps, length(methods), 1),
                  dimnames = list(NULL, methods, "error"))
  
  list(
    settings = list(reps = reps, n = n, n_prob = n_prob, n_nonprob = n_nonprob, K = K,
                    family = family, alpha = alpha, gamma_vec = gamma_vec, pi_vec = pi_vec,
                    Pi = Pi, nonprob_ratio = nonprob_ratio, known_probs = known_probs,
                    val_frac = val_frac, include_simex = include_simex,
                    simex_B = simex_B, simex_method = simex_method),
    true_mean = truth,
    mean_est = mean_est,
    mean_se = mean_se,
    cover = cover,
    errors = errors
  )
}

print_comb_method_summaries <- function(mc_obj) {
  methods <- if (!is.null(mc_obj$mean_est)) colnames(mc_obj$mean_est) else dimnames(mc_obj$errors)[[2]]
  
  cat("\n=== Monte Carlo summary ===\n")
  #print(mc_obj$settings)
  
  if (!is.null(mc_obj$mean_est) && !is.null(mc_obj$mean_se)) {
    out <- t(vapply(methods, function(m) {
      summarize_scalar_mc_with_ci(mc_obj$mean_est[, m], mc_obj$mean_se[, m], mc_obj$true_mean)
    }, numeric(6)))
    out <- data.frame(method = rownames(out), out, row.names = NULL)
    print(out, row.names = FALSE)
    return(invisible(NULL))
  }
  
  errors <- matrix(mc_obj$errors, ncol = length(methods))
  bias <- colMeans(errors, na.rm = TRUE)
  sd   <- apply(errors, 2, sd, na.rm = TRUE)
  mse  <- colMeans(errors^2, na.rm = TRUE)
  
  print(data.frame(method = methods, bias = bias, sd = sd, mse = mse), row.names = F)
  invisible(NULL)
}

# ---------- 5) Convenience: LaTeX tables ----------

latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x
}

format_latex_number <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

format_latex_signif <- function(x, digits = 3) {
  vapply(x, function(z) {
    if (is.na(z)) return("")
    if (z == 0) return(formatC(0, format = "f", digits = digits))
    decimal_places <- if (abs(z) < 1) digits else max(0, digits - floor(log10(abs(z))))
    formatC(z, format = "f", digits = decimal_places)
  }, character(1))
}

latex_param_label <- function(x) {
  x <- as.character(x)
  out <- x
  
  out[x == "gamma"] <- "$\\gamma$"
  out[x == "alpha"] <- "$\\alpha$"
  
  gamma_idx <- grepl("^gamma[0-9]+$", x)
  out[gamma_idx] <- paste0("$\\gamma_{", sub("^gamma", "", x[gamma_idx]), "}$")
  
  alpha_idx <- grepl("^alpha[0-9]+$", x)
  out[alpha_idx] <- paste0("$\\alpha_{", sub("^alpha", "", x[alpha_idx]), "}$")
  
  out
}

format_latex_cell <- function(x, digits = 3, colname = NULL) {
  if (!is.null(colname) && colname == "param") return(latex_param_label(x))
  if (!is.null(colname) && colname == "mse") return(format_latex_signif(1000 * x, digits))
  if (!is.null(colname) && colname %in% c("truth", "bias")) {
    return(ifelse(is.na(x), "", paste0("$", format_latex_number(x, digits), "$")))
  }
  if (is.numeric(x)) return(format_latex_number(x, digits))
  latex_escape(x)
}

latex_header_label <- function(x) {
  if (x == "mse") return("MSE $\\times 10^3$")
  latex_escape(x)
}

format_latex_mass_cell <- function(x, digits = 3, colname = NULL) {
  if (!is.null(colname) && colname == "bias") {
    return(ifelse(is.na(x), "", paste0("$", format_latex_number(1000 * x, digits), "$")))
  }
  if (!is.null(colname) && colname %in% c("sd", "mean_se")) return(format_latex_number(100 * x, digits))
  if (!is.null(colname) && colname == "mse") return(format_latex_signif(10000 * x, digits))
  format_latex_cell(x, digits = digits, colname = colname)
}

latex_mass_header_label <- function(x) {
  if (x == "bias") return("Bias $\\times 10^3$")
  if (x == "sd") return("SD $\\times 10^2$")
  if (x == "mean_se") return("mean\\_se $\\times 10^2$")
  if (x == "mse") return("MSE $\\times 10^4$")
  latex_header_label(x)
}

coefficient_summary_df <- function(mc_obj) {
  est_arr <- mc_obj$est_arr
  true_psi <- mc_obj$true_psi
  methods <- dimnames(est_arr)[[2]]
  pnames <- dimnames(est_arr)[[3]]
  
  out <- list()
  row_id <- 1
  for (m in methods) {
    est_mat <- matrix(est_arr[, m, , drop = FALSE], ncol = length(pnames))
    colnames(est_mat) <- pnames
    
    if (!is.null(mc_obj$se_arr)) {
      se_mat <- matrix(mc_obj$se_arr[, m, , drop = FALSE], ncol = length(pnames))
      colnames(se_mat) <- pnames
      tab <- summarize_mc_with_ci(est_mat, se_mat, true_psi)
    } else {
      tab <- summarize_mc(est_mat, true_psi)
    }
    
    tab <- data.frame(method = m, tab, row.names = NULL)
    out[[row_id]] <- tab
    row_id <- row_id + 1
  }
  
  do.call(rbind, out)
}

mass_summary_df <- function(mc_obj) {
  methods <- if (!is.null(mc_obj$mean_est)) colnames(mc_obj$mean_est) else dimnames(mc_obj$errors)[[2]]
  
  if (!is.null(mc_obj$mean_est) && !is.null(mc_obj$mean_se)) {
    out <- t(vapply(methods, function(m) {
      summarize_scalar_mc_with_ci(mc_obj$mean_est[, m], mc_obj$mean_se[, m], mc_obj$true_mean)
    }, numeric(6)))
    return(data.frame(method = rownames(out), out, row.names = NULL))
  }
  
  errors <- matrix(mc_obj$errors, ncol = length(methods))
  data.frame(
    method = methods,
    bias = colMeans(errors, na.rm = TRUE),
    sd = apply(errors, 2, sd, na.rm = TRUE),
    mse = colMeans(errors^2, na.rm = TRUE),
    row.names = NULL
  )
}

latex_table_from_df <- function(tab, caption = NULL, label = NULL, digits = 3,
                                align = NULL, table_env = TRUE,
                                cell_formatter = format_latex_cell,
                                header_formatter = latex_header_label) {
  if (is.null(align)) {
    align <- paste0("l", paste(rep("r", ncol(tab) - 1), collapse = ""))
  }
  
  formatted <- as.data.frame(setNames(lapply(names(tab), function(nm) {
    cell_formatter(tab[[nm]], digits = digits, colname = nm)
  }), names(tab)), stringsAsFactors = FALSE)
  names(formatted) <- vapply(names(tab), header_formatter, character(1))
  
  lines <- character()
  if (table_env) {
    lines <- c(lines, "\\begin{table}[!ht]", "\\centering")
    if (!is.null(caption)) lines <- c(lines, paste0("\\caption{", latex_escape(caption), "}"))
    if (!is.null(label)) lines <- c(lines, paste0("\\label{", latex_escape(label), "}"))
  }
  
  lines <- c(lines, paste0("\\begin{tabular}{", align, "}"))
  lines <- c(lines, "\\hline")
  lines <- c(lines, paste(names(formatted), collapse = " & "))
  lines <- c(lines, "\\\\")
  lines <- c(lines, "\\hline")
  
  for (i in seq_len(nrow(formatted))) {
    lines <- c(lines, paste(as.character(formatted[i, ]), collapse = " & "))
    lines <- c(lines, "\\\\")
  }
  
  lines <- c(lines, "\\hline", "\\end{tabular}")
  if (table_env) lines <- c(lines, "\\end{table}")
  
  paste(lines, collapse = "\n")
}

print_latex_method_summaries <- function(mc_obj, digits = 3, caption = NULL,
                                         label = NULL, file = NULL,
                                         table_env = TRUE) {
  tex <- latex_table_from_df(coefficient_summary_df(mc_obj), caption = caption,
                             label = label, digits = digits, table_env = table_env)
  if (is.null(file)) {
    cat(tex, "\n")
  } else {
    writeLines(tex, con = file)
  }
  invisible(tex)
}

print_latex_comb_method_summaries <- function(mc_obj, digits = 3, caption = NULL,
                                              label = NULL, file = NULL,
                                              table_env = TRUE) {
  tex <- latex_table_from_df(mass_summary_df(mc_obj), caption = caption,
                             label = label, digits = digits, table_env = table_env,
                             cell_formatter = format_latex_mass_cell,
                             header_formatter = latex_mass_header_label)
  if (is.null(file)) {
    cat(tex, "\n")
  } else {
    writeLines(tex, con = file)
  }
  invisible(tex)
}

print_latex_summaries <- function(mc_obj, digits = 3, caption = NULL,
                                  label = NULL, file = NULL,
                                  table_env = TRUE) {
  if (!is.null(mc_obj$est_arr)) {
    return(print_latex_method_summaries(mc_obj, digits = digits, caption = caption,
                                        label = label, file = file,
                                        table_env = table_env))
  }
  if (!is.null(mc_obj$errors) || !is.null(mc_obj$mean_est)) {
    return(print_latex_comb_method_summaries(mc_obj, digits = digits, caption = caption,
                                             label = label, file = file,
                                             table_env = table_env))
  }
  stop("Unsupported simulation object.")
}


# ============================================================
# Simulations:
mc1 <- run_comb_binary(reps = 1000, n = 100000, n_nonprob = 10000, n_prob = 2000,
                       p01 = 0.15, p10 = 0.1, known_probs = FALSE, val_frac = 0.1)
mc2 <- run_comb_binary(reps = 1000, n = 100000, n_nonprob = 10000, n_prob = 2000,
                       p01 = 0.2, p10 = 0.15, known_probs = FALSE, val_frac = 0.02)
mc3 <- run_comb_binary(reps = 1000, n = 100000, n_nonprob = 10000, n_prob = 2000, gamma = 1.3, alpha = c(-1, 1.5),
                       p01 = 0.4, p10 = 0.2, known_probs = FALSE, val_frac = 0.1)
print_comb_method_summaries(mc1)
print_comb_method_summaries(mc2)
print_comb_method_summaries(mc3)
# print_latex_summaries(mc1)
# print_latex_summaries(mc2)
# print_latex_summaries(mc3)
