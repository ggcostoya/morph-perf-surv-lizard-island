# Priors for the model in 02-model.qmd. Shared by the simulator and the
# Stan model so the two cannot diverge.
#
# Hierarchical: species values ~ normal(mu, tau), mu ~ normal(mu_mean,
# mu_sd), tau ~ half-normal(0, tau_sd); one mu and tau per trait.
priors <- list(
  # performance
  alpha = list(sd = 0.5),
  delta = list(mu_mean = 0, mu_sd = 1, tau_sd = 0.5),
  b_svl = list(mu_mean = 0, mu_sd = 0.5, tau_sd = 0.25),
  b_shape = list(mu_mean = 0, mu_sd = 0.5, tau_sd = 0.25),
  t = list(mu_mean = 0, mu_sd = 0.5, tau_sd = 0.25),
  sigma = list(sd = 1),
  Omega = list(eta = 2),
  # survival
  a = list(mu_mean = -1, mu_sd = 1, tau_sd = 0.5),
  c = list(mu_mean = 0, mu_sd = 1, tau_sd = 0.5),
  beta = list(mu_mean = 0, mu_sd = 0.5, tau_sd = 0.25),
  q = list(mu_mean = 0, mu_sd = 0.25, tau_sd = 0.1),
  gamma = list(mu_mean = 0, mu_sd = 0.5, tau_sd = 0.25),
  # detection
  pi = list(mu_mean = 0, mu_sd = 1.5, tau_sd = 1),
  ell = list(mu_mean = 0, mu_sd = 1.5, tau_sd = 1),
  zeta = list(sd = 1)
)
