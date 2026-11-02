# Assembles the data list for stan/model.stan from a design or data list
# (prepare_data() or simulate_data()) and the priors (R/priors.R).

# priors$delta$mu_sd -> delta_mu_sd, priors$alpha$sd -> alpha_sd, ...
flatten_priors <- function(priors) {
  unlist(lapply(names(priors), \(p) {
    setNames(as.list(unlist(priors[[p]])), paste(p, names(priors[[p]]), sep = "_"))
  }), recursive = FALSE)
}

make_stan_data <- function(data, priors, direct_effect = FALSE) {
  c(
    data[c(
      "N", "S", "J", "K", "species", "male", "x", "w", "shape",
      "u", "obs", "cens", "c_leaf", "y2", "y3"
    )],
    list(
      J_t = ncol(data$w),
      temp_idx = match(temp_traits, perf_traits),
      leaf = match("leaf", perf_traits),
      direct_effect = as.integer(direct_effect)
    ),
    flatten_priors(priors)
  )
}
