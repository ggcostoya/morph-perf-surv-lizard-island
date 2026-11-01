# Simulator for the model in 02-model.qmd. Conditions on the observed
# design from prepare_data() (species, sex, morphology, trial temperatures,
# missing-performance pattern) and simulates performance, survival and
# resighting from parameters drawn from the priors (R/priors.R).

## Helpers ----

# S x n matrix of species values with one mu and tau per column
draw_hier <- function(S, n, prior) {
  mu <- rnorm(n, prior$mu_mean, prior$mu_sd)
  tau <- abs(rnorm(n, 0, prior$tau_sd))
  value <- matrix(rnorm(S * n, rep(mu, each = S), rep(tau, each = S)), S, n)
  list(mu = mu, tau = tau, value = value)
}

# LKJ correlation matrix draw (onion method)
rlkjcorr <- function(d, eta) {
  if (d == 1) return(matrix(1))
  alpha <- eta + (d - 2) / 2
  r12 <- 2 * rbeta(1, alpha, alpha) - 1
  R <- matrix(0, d, d)
  R[1, 1] <- 1
  R[1, 2] <- r12
  R[2, 2] <- sqrt(1 - r12^2)
  if (d > 2) {
    for (m in 2:(d - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m / 2, alpha)
      z <- rnorm(m)
      z <- z / sqrt(sum(z^2))
      R[1:m, m + 1] <- sqrt(y) * z
      R[m + 1, m + 1] <- sqrt(1 - y)
    }
  }
  crossprod(R)
}

# flatten a parameter into Stan-style names, e.g. beta[3,2]
stan_flatten <- function(x, name) {
  if (length(x) == 1 && is.null(dim(x))) return(setNames(x, name))
  idx <- arrayInd(seq_along(x), dim(x) %||% length(x))
  setNames(as.vector(x), paste0(name, "[", apply(idx, 1, paste, collapse = ","), "]"))
}

## Draw parameters from the prior ----

draw_params <- function(design, priors, direct_effect = FALSE) {
  S <- design$S
  J <- design$J
  K <- design$K
  J_t <- ncol(design$w)

  delta <- draw_hier(S, J, priors$delta)
  b_svl <- draw_hier(S, J, priors$b_svl)
  b_shape <- draw_hier(S, J, priors$b_shape)
  t <- draw_hier(S, J_t, priors$t)
  a <- draw_hier(S, 1, priors$a)
  c <- draw_hier(S, 1, priors$c)
  beta <- draw_hier(S, J, priors$beta)
  q <- draw_hier(S, J, priors$q)
  pi <- draw_hier(S, 1, priors$pi)
  ell <- draw_hier(S, 1, priors$ell)

  params <- list(
    # performance
    alpha = matrix(rnorm(S * J, 0, priors$alpha$sd), S, J),
    mu_delta = delta$mu, tau_delta = delta$tau, delta = delta$value,
    mu_b_svl = b_svl$mu, tau_b_svl = b_svl$tau, b_svl = b_svl$value,
    mu_b_shape = b_shape$mu, tau_b_shape = b_shape$tau, b_shape = b_shape$value,
    mu_t = t$mu, tau_t = t$tau, t = t$value,
    sigma = matrix(abs(rnorm(S * J, 0, priors$sigma$sd)), S, J),
    Omega = rlkjcorr(J, priors$Omega$eta),
    # survival
    mu_a = a$mu, tau_a = a$tau, a = as.vector(a$value),
    mu_c = c$mu, tau_c = c$tau, c = as.vector(c$value),
    mu_beta = beta$mu, tau_beta = beta$tau, beta = beta$value,
    mu_q = q$mu, tau_q = q$tau, q = q$value,
    # detection
    mu_pi = pi$mu, tau_pi = pi$tau, pi = as.vector(pi$value),
    mu_ell = ell$mu, tau_ell = ell$tau, ell = as.vector(ell$value),
    zeta_p = rnorm(1, 0, priors$zeta$sd),
    zeta_lambda = rnorm(1, 0, priors$zeta$sd)
  )

  # direct morphology -> survival; absent (0) under full mediation (M1)
  if (direct_effect) {
    gamma <- draw_hier(S, K, priors$gamma)
    params$mu_gamma <- gamma$mu
    params$tau_gamma <- gamma$tau
    params$gamma <- gamma$value
  } else {
    params$gamma <- matrix(0, S, K)
  }
  params
}

## Simulate outcomes given parameters and the observed design ----

simulate_data <- function(params, design) {
  N <- design$N
  J <- design$J
  s <- design$species
  m <- design$male
  x <- design$x
  temp_cols <- match(temp_traits, perf_traits)

  # performance: u_i ~ MVN(mu_i, Sigma_s)
  temp <- matrix(0, N, J)
  temp[, temp_cols] <- params$t[s, , drop = FALSE] * design$w
  mu <- params$alpha[s, , drop = FALSE] +
    params$delta[s, , drop = FALSE] * m +
    params$b_svl[s, , drop = FALSE] * x[, "svl"] +
    params$b_shape[s, , drop = FALSE] * x[, design$shape] +
    temp
  u <- matrix(NA_real_, N, J, dimnames = list(NULL, perf_traits))
  for (sp in seq_len(design$S)) {
    rows <- s == sp
    if (!any(rows)) next
    sd_sp <- diag(params$sigma[sp, ], J)
    u[rows, ] <- mu[rows, , drop = FALSE] +
      MASS::mvrnorm(sum(rows), rep(0, J), sd_sp %*% params$Omega %*% sd_sp)
  }

  # survival depends on temperature-adjusted performance
  u_adj <- u - temp
  eta <- params$a[s] + params$c[s] * m +
    rowSums(u_adj * params$beta[s, , drop = FALSE]) +
    rowSums(u_adj^2 * params$q[s, , drop = FALSE]) +
    rowSums(x * params$gamma[s, , drop = FALSE])
  z <- rbinom(N, 1, plogis(eta))

  # resighting in fall 2024 (p) and spring 2025 (lambda)
  p <- plogis(params$pi[s] + params$zeta_p * m)
  lambda <- plogis(params$ell[s] + params$zeta_lambda * m)
  y2 <- rbinom(N, 1, z * p)
  y3 <- rbinom(N, 1, z * lambda)

  # observation: copy the real missing pattern, censor leaf below limit
  leaf <- match("leaf", perf_traits)
  cens <- as.integer(design$obs[, leaf] == 1 & u[, leaf] < design$c_leaf[s])
  u_obs <- u
  u_obs[design$obs == 0] <- 0
  u_obs[cens == 1, leaf] <- 0

  modifyList(design, list(u = u_obs, cens = cens, y2 = y2, y3 = y3))
}

## One SBC draw: true parameters + model-ready data ----

# Matches SBC::SBC_generator_function(); wrap variables in
# posterior::as_draws_matrix(t(...)) when running with SBC.
simulate_one <- function(design, priors, direct_effect = FALSE) {
  params <- draw_params(design, priors, direct_effect)
  generated <- simulate_data(params, design)

  # report only parameters the model estimates
  report <- params
  report$Omega <- NULL
  if (!direct_effect) report$gamma <- NULL
  variables <- unlist(unname(imap(report, stan_flatten)))
  omega <- stan_flatten(params$Omega, "Omega")
  variables <- c(variables, omega[lower.tri(params$Omega)])

  list(variables = variables, generated = generated)
}
