
## Packages ----

library(tidyverse)

## Study design ----

# Morphology is observed and conditioned on, so its distribution is part of
# the design rather than the prior. All values are placeholders.
design <- list(
  species = c("caro", "cris", "sag"),
  n_individuals = c(80, 500, 80),
  morph_traits = c("svl", "hindlimb"),
  perf_traits = c("sprint", "cling"),
  # species x morphology trait
  morph_means = rbind(
    caro = c(svl = 5, hindlimb = 3),
    cris = c(svl = 10, hindlimb = 6),
    sag = c(svl = 9, hindlimb = 5)
  ),
  morph_sds = rbind(
    caro = c(svl = 1, hindlimb = 1),
    cris = c(svl = 1, hindlimb = 1),
    sag = c(svl = 1, hindlimb = 1)
  ),
  # correlation among morphology traits (shared across species)
  morph_cor = matrix(c(1, 0.7, 0.7, 1), 2, 2)
)

## Priors ----

# Hierarchical parameters: species values ~ normal(mu, tau), with
# mu ~ normal(mu_mean, mu_sd) and tau ~ half-normal(0, tau_sd).
# Each entry can be a scalar or one value per trait (or trait pair, in
# column-major order). The Stan model must use exactly these priors.
priors <- list(
  alpha_p = list(mu_mean = 4, mu_sd = 1, tau_sd = 0.5),     # perf at mean morph
  B_p = list(mu_mean = 0.4, mu_sd = 0.2, tau_sd = 0.1),     # perf ~ morph slopes
  sigma_p = list(sd = 0.5),                                 # perf residual SD
  Omega_p = list(eta = 2),                                  # perf residual correlation (LKJ)
  alpha_s = list(mu_mean = -1.5, mu_sd = 0.5, tau_sd = 0.5), # logit surv at mean perf
  beta_s = list(mu_mean = 0.6, mu_sd = 0.3, tau_sd = 0.2),   # linear selection on perf
  q_s = list(mu_mean = 0, mu_sd = 0.2, tau_sd = 0.1),        # quadratic selection on perf
  gamma_s = list(mu_mean = 0, mu_sd = 0.3, tau_sd = 0.1)     # direct morph -> surv
)

## Helpers ----

# species-level array (dims = c(S, ...)) with one mu and tau per trait element
draw_hier <- function(dims, prior) {
  n_el <- prod(dims[-1])
  mu <- rnorm(n_el, prior$mu_mean, prior$mu_sd)
  tau <- abs(rnorm(n_el, 0, prior$tau_sd))
  value <- matrix(
    rnorm(dims[1] * n_el, rep(mu, each = dims[1]), rep(tau, each = dims[1])),
    dims[1], n_el
  )
  if (length(dims) > 1) {
    mu <- array(mu, dims[-1])
    tau <- array(tau, dims[-1])
  }
  list(mu = mu, tau = tau, value = array(value, dims))
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

# subtract species-specific column means
center_by <- function(x, s) {
  for (sp in unique(s)) {
    rows <- s == sp
    x[rows, ] <- sweep(x[rows, , drop = FALSE], 2, colMeans(x[rows, , drop = FALSE]))
  }
  x
}

# flatten a parameter into Stan-style names, e.g. B_p[3,1,2]
stan_flatten <- function(x, name) {
  if (length(x) == 1 && is.null(dim(x))) return(setNames(x, name))
  idx <- arrayInd(seq_along(x), dim(x) %||% length(x))
  setNames(as.vector(x), paste0(name, "[", apply(idx, 1, paste, collapse = ","), "]"))
}

## Draw parameters from the prior ----

draw_params <- function(design, priors, direct_effect = FALSE) {
  S <- length(design$species)
  K <- length(design$morph_traits)
  J <- length(design$perf_traits)

  alpha_p <- draw_hier(c(S, J), priors$alpha_p)
  B_p <- draw_hier(c(S, J, K), priors$B_p)
  alpha_s <- draw_hier(S, priors$alpha_s)
  beta_s <- draw_hier(c(S, J), priors$beta_s)
  q_s <- draw_hier(c(S, J), priors$q_s)

  params <- list(
    mu_alpha_p = alpha_p$mu, tau_alpha_p = alpha_p$tau, alpha_p = alpha_p$value,
    mu_B_p = B_p$mu, tau_B_p = B_p$tau, B_p = B_p$value,
    sigma_p = matrix(abs(rnorm(S * J, 0, priors$sigma_p$sd)), S, J),
    Omega_p = rlkjcorr(J, priors$Omega_p$eta),
    mu_alpha_s = alpha_s$mu, tau_alpha_s = alpha_s$tau, alpha_s = alpha_s$value,
    mu_beta_s = beta_s$mu, tau_beta_s = beta_s$tau, beta_s = beta_s$value,
    mu_q_s = q_s$mu, tau_q_s = q_s$tau, q_s = q_s$value
  )

  # direct morphology -> survival path; absent (0) under full mediation
  if (direct_effect) {
    gamma_s <- draw_hier(c(S, K), priors$gamma_s)
    params$mu_gamma_s <- gamma_s$mu
    params$tau_gamma_s <- gamma_s$tau
    params$gamma_s <- gamma_s$value
  } else {
    params$gamma_s <- matrix(0, S, K)
  }
  params
}

## Simulate one data set given parameters ----

# Morphology and performance are centered on species sample means, so
# alpha_p is mean performance and plogis(alpha_s) is survival of an
# average individual.
simulate_data <- function(params, design) {
  S <- length(design$species)
  K <- length(design$morph_traits)
  J <- length(design$perf_traits)
  s <- rep(seq_len(S), design$n_individuals)
  N <- length(s)

  # morphology
  morph <- matrix(NA_real_, N, K)
  for (sp in seq_len(S)) {
    sd_m <- diag(design$morph_sds[sp, ], K)
    morph[s == sp, ] <- MASS::mvrnorm(
      design$n_individuals[sp], design$morph_means[sp, ], sd_m %*% design$morph_cor %*% sd_m
    )
  }
  morph_c <- center_by(morph, s)

  # performance (measured without error; residual = biological variation)
  perf <- matrix(NA_real_, N, J)
  for (sp in seq_len(S)) {
    rows <- s == sp
    B <- matrix(params$B_p[sp, , ], J, K)
    sd_p <- diag(params$sigma_p[sp, ], J)
    resid <- MASS::mvrnorm(design$n_individuals[sp], rep(0, J), sd_p %*% params$Omega_p %*% sd_p)
    perf[rows, ] <- sweep(morph_c[rows, , drop = FALSE] %*% t(B), 2, params$alpha_p[sp, ], "+") +
      matrix(resid, ncol = J)
  }
  perf_c <- center_by(perf, s)

  # survival between the two censuses
  eta <- params$alpha_s[s] +
    rowSums(perf_c * params$beta_s[s, , drop = FALSE]) +
    rowSums(perf_c^2 * params$q_s[s, , drop = FALSE]) +
    rowSums(morph_c * params$gamma_s[s, , drop = FALSE])
  surv <- rbinom(N, 1, plogis(eta))

  colnames(morph) <- colnames(morph_c) <- design$morph_traits
  colnames(perf) <- colnames(perf_c) <- design$perf_traits
  list(s = s, morph = morph, morph_c = morph_c, perf = perf, perf_c = perf_c, surv = surv)
}

# tidy data frame for plotting
as_tibble_sim <- function(sim, design) {
  tibble(species = factor(design$species[sim$s], levels = design$species)) |>
    bind_cols(as_tibble(sim$morph), as_tibble(sim$perf)) |>
    mutate(surv = sim$surv)
}

## One SBC draw: true parameters + model-ready data ----

# Matches SBC::SBC_generator_function(); wrap variables in
# posterior::as_draws_matrix(t(...)) when running with SBC.
simulate_one <- function(design, priors, direct_effect = FALSE) {
  params <- draw_params(design, priors, direct_effect)
  sim <- simulate_data(params, design)

  # report only parameters the model estimates; drop the fixed
  # diagonal of Omega_p and gamma_s under full mediation
  J <- length(design$perf_traits)
  report <- params
  report$Omega_p <- NULL
  if (!direct_effect) report$gamma_s <- NULL
  variables <- unlist(unname(imap(report, stan_flatten)))
  if (J > 1) {
    omega <- stan_flatten(params$Omega_p, "Omega_p")
    variables <- c(variables, omega[lower.tri(params$Omega_p)])
  }

  list(
    variables = variables,
    generated = list(
      N = length(sim$s),
      S = length(design$species),
      K = length(design$morph_traits),
      J = J,
      species = sim$s,
      morph_c = sim$morph_c,
      perf = sim$perf,
      perf_c = sim$perf_c,
      surv = sim$surv
    )
  )
}

## Test ----

set.seed(1)
params <- draw_params(design, priors, direct_effect = TRUE)
test_data <- simulate_data(params, design) |> as_tibble_sim(design)

test_data |>
  pivot_longer(all_of(design$morph_traits), names_to = "morph_trait", values_to = "morph") |>
  pivot_longer(all_of(design$perf_traits), names_to = "perf_trait", values_to = "perf") |>
  ggplot(aes(x = morph, y = perf, color = species)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm") +
  facet_grid(perf_trait ~ morph_trait, scales = "free")

test_data |>
  pivot_longer(all_of(design$perf_traits), names_to = "perf_trait", values_to = "perf") |>
  ggplot(aes(x = perf, y = surv, color = species)) +
  geom_point(alpha = 0.3) +
  geom_smooth(
    method = "glm",
    formula = y ~ poly(x, 2),
    method.args = list(family = "binomial")
  ) +
  facet_wrap(~perf_trait, scales = "free_x")

# prior predictive check: survival rate per species across draws
prior_surv <- map(1:200, \(i) {
  d <- simulate_one(design, priors)$generated
  tapply(d$surv, d$species, mean)
}) |>
  do.call(what = rbind)
colnames(prior_surv) <- design$species
summary(prior_surv)
