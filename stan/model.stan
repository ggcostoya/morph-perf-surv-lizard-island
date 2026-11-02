// Morphology -> performance -> survival, with imperfect detection.
// Specification: 02-model.qmd. Priors are passed as data from R/priors.R.
// direct_effect = 0 gives M1 (full mediation), 1 gives M2 (partial).

data {
  int<lower=1> N; // lizards
  int<lower=1> S; // species
  int<lower=1> J; // performance traits
  int<lower=1> K; // morphological predictors (column 1 = size)
  int<lower=0> J_t; // performance traits with a temperature term
  array[N] int<lower=1, upper=S> species;
  vector<lower=0, upper=1>[N] male;
  matrix[N, K] x; // standardised morphology
  matrix[N, J_t] w; // standardised trial temperature
  array[J] int<lower=2, upper=K> shape; // shape trait k(j) for each j
  array[J_t] int<lower=1, upper=J> temp_idx; // traits with temperature
  int<lower=1, upper=J> leaf; // trait with left-censoring

  matrix[N, J] u; // standardised log performance (placeholders where unobserved)
  array[N, J] int<lower=0, upper=1> obs; // measured (incl. censored)
  array[N] int<lower=0, upper=1> cens; // leaf value below detection limit
  vector[S] c_leaf; // detection limit, standardised scale

  array[N] int<lower=0, upper=1> y2; // seen fall 2024
  array[N] int<lower=0, upper=1> y3; // seen spring 2025

  int<lower=0, upper=1> direct_effect;

  // priors
  real<lower=0> alpha_sd;
  real delta_mu_mean;
  real<lower=0> delta_mu_sd;
  real<lower=0> delta_tau_sd;
  real b_svl_mu_mean;
  real<lower=0> b_svl_mu_sd;
  real<lower=0> b_svl_tau_sd;
  real b_shape_mu_mean;
  real<lower=0> b_shape_mu_sd;
  real<lower=0> b_shape_tau_sd;
  real t_mu_mean;
  real<lower=0> t_mu_sd;
  real<lower=0> t_tau_sd;
  real log_t2_mu_mean;
  real<lower=0> log_t2_mu_sd;
  real<lower=0> log_t2_tau_sd;
  real<lower=0> sigma_sd;
  real<lower=0> Omega_eta;
  real a_mu_mean;
  real<lower=0> a_mu_sd;
  real<lower=0> a_tau_sd;
  real c_mu_mean;
  real<lower=0> c_mu_sd;
  real<lower=0> c_tau_sd;
  real beta_mu_mean;
  real<lower=0> beta_mu_sd;
  real<lower=0> beta_tau_sd;
  real q_mu_mean;
  real<lower=0> q_mu_sd;
  real<lower=0> q_tau_sd;
  real gamma_mu_mean;
  real<lower=0> gamma_mu_sd;
  real<lower=0> gamma_tau_sd;
  real pi_mu_mean;
  real<lower=0> pi_mu_sd;
  real<lower=0> pi_tau_sd;
  real ell_mu_mean;
  real<lower=0> ell_mu_sd;
  real<lower=0> ell_tau_sd;
  real<lower=0> zeta_sd;
}

transformed data {
  // positions of missing and censored performance values
  int N_mis = 0;
  for (i in 1:N) {
    for (j in 1:J) {
      N_mis += 1 - obs[i, j];
    }
  }
  int N_cens = sum(cens);
  array[N_mis] int mis_i;
  array[N_mis] int mis_j;
  array[N_cens] int cens_i;
  vector[N_cens] cens_upper;
  {
    int m = 1;
    int c = 1;
    for (i in 1:N) {
      for (j in 1:J) {
        if (obs[i, j] == 0) {
          mis_i[m] = i;
          mis_j[m] = j;
          m += 1;
        }
      }
      if (cens[i] == 1) {
        cens_i[c] = i;
        cens_upper[c] = c_leaf[species[i]];
        c += 1;
      }
    }
  }
  // lizards grouped by species (for the performance likelihood)
  array[S] int n_sp = rep_array(0, S);
  for (i in 1:N) {
    n_sp[species[i]] += 1;
  }
  array[N] int by_sp;
  array[S] int sp_start;
  {
    array[S] int pos;
    sp_start[1] = 1;
    for (s in 2:S) {
      sp_start[s] = sp_start[s - 1] + n_sp[s - 1];
    }
    pos = sp_start;
    for (i in 1:N) {
      by_sp[pos[species[i]]] = i;
      pos[species[i]] += 1;
    }
  }

  // lizards grouped by resighting history (y2 y3)
  int n_11 = 0;
  int n_01 = 0;
  int n_10 = 0;
  for (i in 1:N) {
    n_11 += y2[i] * y3[i];
    n_01 += (1 - y2[i]) * y3[i];
    n_10 += y2[i] * (1 - y3[i]);
  }
  int n_00 = N - n_11 - n_01 - n_10;
  array[n_11] int idx_11;
  array[n_01] int idx_01;
  array[n_10] int idx_10;
  array[n_00] int idx_00;
  {
    array[4] int k = {1, 1, 1, 1};
    for (i in 1:N) {
      if (y2[i] == 1 && y3[i] == 1) {
        idx_11[k[1]] = i;
        k[1] += 1;
      } else if (y3[i] == 1) {
        idx_01[k[2]] = i;
        k[2] += 1;
      } else if (y2[i] == 1) {
        idx_10[k[3]] = i;
        k[3] += 1;
      } else {
        idx_00[k[4]] = i;
        k[4] += 1;
      }
    }
  }
  matrix[N, J] male_J = rep_matrix(male, J);
  matrix[N, J] svl_J = rep_matrix(x[, 1], J);
  matrix[N, J] x_shape = x[, shape];
  int K_g = direct_effect ? K : 0;
  int S_g = direct_effect ? S : 0;
}

parameters {
  // performance
  matrix[S, J] alpha;
  vector[J] mu_delta;
  vector<lower=0>[J] tau_delta;
  matrix[S, J] z_delta;
  vector[J] mu_b_svl;
  vector<lower=0>[J] tau_b_svl;
  matrix[S, J] z_b_svl;
  vector[J] mu_b_shape;
  vector<lower=0>[J] tau_b_shape;
  matrix[S, J] z_b_shape;
  vector[J_t] mu_t;
  vector<lower=0>[J_t] tau_t;
  matrix[S, J_t] z_t;
  vector[J_t] mu_log_t2; // log curvature of the temperature effect
  vector<lower=0>[J_t] tau_log_t2;
  matrix[S, J_t] z_t2;
  matrix<lower=0>[S, J] sigma;
  cholesky_factor_corr[J] L_Omega;
  vector[N_mis] u_mis;
  vector<upper=cens_upper>[N_cens] u_cens;

  // survival
  real mu_a;
  real<lower=0> tau_a;
  vector[S] z_a;
  real mu_c;
  real<lower=0> tau_c;
  vector[S] z_c;
  vector[J] mu_beta;
  vector<lower=0>[J] tau_beta;
  matrix[S, J] z_beta;
  vector[J] mu_q;
  vector<lower=0>[J] tau_q;
  matrix[S, J] z_q;
  vector[K_g] mu_gamma;
  vector<lower=0>[K_g] tau_gamma;
  matrix[S_g, K_g] z_gamma;

  // detection
  real mu_pi;
  real<lower=0> tau_pi;
  vector[S] z_pi;
  real mu_ell;
  real<lower=0> tau_ell;
  vector[S] z_ell;
  real zeta_p;
  real zeta_lambda;
}

transformed parameters {
  // non-centred species-level parameters
  matrix[S, J] delta = rep_matrix(mu_delta', S) + diag_post_multiply(z_delta, tau_delta);
  matrix[S, J] b_svl = rep_matrix(mu_b_svl', S) + diag_post_multiply(z_b_svl, tau_b_svl);
  matrix[S, J] b_shape = rep_matrix(mu_b_shape', S) + diag_post_multiply(z_b_shape, tau_b_shape);
  matrix[S, J_t] t = rep_matrix(mu_t', S) + diag_post_multiply(z_t, tau_t);
  // quadratic temperature effect, concave by construction
  matrix[S, J_t] t2 = -exp(rep_matrix(mu_log_t2', S) + diag_post_multiply(z_t2, tau_log_t2));
  vector[S] a = mu_a + tau_a * z_a;
  vector[S] c = mu_c + tau_c * z_c;
  matrix[S, J] beta = rep_matrix(mu_beta', S) + diag_post_multiply(z_beta, tau_beta);
  matrix[S, J] q = rep_matrix(mu_q', S) + diag_post_multiply(z_q, tau_q);
  matrix[S_g, K_g] gamma = rep_matrix(mu_gamma', S_g) + diag_post_multiply(z_gamma, tau_gamma);
  vector[S] pi = mu_pi + tau_pi * z_pi;
  vector[S] ell = mu_ell + tau_ell * z_ell;
}

model {
  // complete performance matrix: observed, missing and censored values
  matrix[N, J] u_full = u;
  for (m in 1:N_mis) {
    u_full[mis_i[m], mis_j[m]] = u_mis[m];
  }
  for (n in 1:N_cens) {
    u_full[cens_i[n], leaf] = u_cens[n];
  }

  // temperature contribution to each performance value (quadratic)
  matrix[N, J] temp = rep_matrix(0, N, J);
  for (jt in 1:J_t) {
    temp[, temp_idx[jt]] = to_vector(t[species, jt]) .* w[, jt]
                           + to_vector(t2[species, jt]) .* square(w[, jt]);
  }

  // performance submodel: u_i ~ MVN(mu_i, Sigma_s), one solve per species
  matrix[N, J] mu = alpha[species] + delta[species] .* male_J
                    + b_svl[species] .* svl_J + b_shape[species] .* x_shape
                    + temp;
  for (s in 1:S) {
    array[n_sp[s]] int sp_idx = by_sp[sp_start[s]:(sp_start[s] + n_sp[s] - 1)];
    matrix[J, J] L_Sigma = diag_pre_multiply(sigma[s]', L_Omega);
    matrix[J, n_sp[s]] z = mdivide_left_tri_low(L_Sigma, (u_full[sp_idx] - mu[sp_idx])');
    target += std_normal_lpdf(to_vector(z)) - n_sp[s] * sum(log(diagonal(L_Sigma)));
  }

  // survival on temperature-adjusted performance
  matrix[N, J] u_adj = u_full - temp;
  vector[N] eta = a[species] + c[species] .* male
                  + rows_dot_product(u_adj, beta[species])
                  + rows_dot_product(square(u_adj), q[species]);
  if (direct_effect) {
    eta += rows_dot_product(x, gamma[species]);
  }

  // detection: fall (p) and spring 2025 (lambda); likelihood by history
  vector[N] p = inv_logit(pi[species] + zeta_p * male);
  vector[N] lambda = inv_logit(ell[species] + zeta_lambda * male);
  target += sum(log_inv_logit(eta[idx_11]) + log(p[idx_11]) + log(lambda[idx_11]));
  target += sum(log_inv_logit(eta[idx_01]) + log1m(p[idx_01]) + log(lambda[idx_01]));
  target += sum(log_inv_logit(eta[idx_10]) + log(p[idx_10]) + log1m(lambda[idx_10]));
  // 00: dead, or alive and missed twice
  target += sum(log1m(inv_logit(eta[idx_00])
                      .* (1 - (1 - p[idx_00]) .* (1 - lambda[idx_00]))));

  // priors: performance
  to_vector(alpha) ~ normal(0, alpha_sd);
  mu_delta ~ normal(delta_mu_mean, delta_mu_sd);
  tau_delta ~ normal(0, delta_tau_sd);
  to_vector(z_delta) ~ std_normal();
  mu_b_svl ~ normal(b_svl_mu_mean, b_svl_mu_sd);
  tau_b_svl ~ normal(0, b_svl_tau_sd);
  to_vector(z_b_svl) ~ std_normal();
  mu_b_shape ~ normal(b_shape_mu_mean, b_shape_mu_sd);
  tau_b_shape ~ normal(0, b_shape_tau_sd);
  to_vector(z_b_shape) ~ std_normal();
  mu_t ~ normal(t_mu_mean, t_mu_sd);
  tau_t ~ normal(0, t_tau_sd);
  to_vector(z_t) ~ std_normal();
  mu_log_t2 ~ normal(log_t2_mu_mean, log_t2_mu_sd);
  tau_log_t2 ~ normal(0, log_t2_tau_sd);
  to_vector(z_t2) ~ std_normal();
  to_vector(sigma) ~ normal(0, sigma_sd);
  L_Omega ~ lkj_corr_cholesky(Omega_eta);

  // priors: survival
  mu_a ~ normal(a_mu_mean, a_mu_sd);
  tau_a ~ normal(0, a_tau_sd);
  z_a ~ std_normal();
  mu_c ~ normal(c_mu_mean, c_mu_sd);
  tau_c ~ normal(0, c_tau_sd);
  z_c ~ std_normal();
  mu_beta ~ normal(beta_mu_mean, beta_mu_sd);
  tau_beta ~ normal(0, beta_tau_sd);
  to_vector(z_beta) ~ std_normal();
  mu_q ~ normal(q_mu_mean, q_mu_sd);
  tau_q ~ normal(0, q_tau_sd);
  to_vector(z_q) ~ std_normal();
  mu_gamma ~ normal(gamma_mu_mean, gamma_mu_sd);
  tau_gamma ~ normal(0, gamma_tau_sd);
  to_vector(z_gamma) ~ std_normal();

  // priors: detection
  mu_pi ~ normal(pi_mu_mean, pi_mu_sd);
  tau_pi ~ normal(0, pi_tau_sd);
  z_pi ~ std_normal();
  mu_ell ~ normal(ell_mu_mean, ell_mu_sd);
  tau_ell ~ normal(0, ell_tau_sd);
  z_ell ~ std_normal();
  zeta_p ~ normal(0, zeta_sd);
  zeta_lambda ~ normal(0, zeta_sd);
}

generated quantities {
  corr_matrix[J] Omega = multiply_lower_tri_self_transpose(L_Omega);

  // effects of each morphological predictor on logit survival at the
  // species mean: indirect (through performance), direct and total
  matrix[S, K] indirect = rep_matrix(0, S, K);
  for (s in 1:S) {
    for (j in 1:J) {
      indirect[s, 1] += b_svl[s, j] * beta[s, j];
      indirect[s, shape[j]] += b_shape[s, j] * beta[s, j];
    }
  }
  matrix[S, K] total = indirect;
  if (direct_effect) {
    total += gamma;
  }
}
