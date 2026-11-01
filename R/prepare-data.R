# Turns data/lizards-spring-2024.csv into the model's data list. Used for
# the real fit and, through its design columns, by the simulator.

species_levels <- c("caro", "cris", "sag")
morph_traits <- c("svl", "hl", "hw", "toe", "pad")
perf_traits <- c("sprint", "bite", "bark", "leaf")
# morphological shape trait linked to each performance trait, k(j)
shape_of <- c(sprint = "hl", bite = "hw", bark = "toe", leaf = "pad")
# performance traits with a body-temperature term
temp_traits <- c("sprint", "bite")
leaf_limit <- 0.01

# standardise within species using the observed (non-missing) values
std_by <- function(v, s) {
  for (sp in unique(s)) {
    rows <- s == sp
    v[rows] <- (v[rows] - mean(v[rows], na.rm = TRUE)) / sd(v[rows], na.rm = TRUE)
  }
  v
}

prepare_data <- function(data) {
  d <- data |>
    dplyr::filter(is.na(exclusion)) |>
    dplyr::mutate(species = factor(species, levels = species_levels))
  s <- as.integer(d$species)

  lsvl <- log(d$svl)
  x <- cbind(
    svl = std_by(lsvl, s),
    hl = std_by(log(d$hindlimb) - lsvl, s),
    hw = std_by(log(d$head_width) - lsvl, s),
    toe = std_by(log(d$hind_toe) - lsvl, s),
    pad = std_by(log(d$toepad_area) - 2 * lsvl, s)
  )

  # trial temperatures; missing set to the species mean (0)
  w <- cbind(sprint = std_by(d$tb_sprint, s), bite = std_by(d$tb_bite, s))
  w[is.na(w)] <- 0

  # leaf zeros are censored below the limit; standardise on non-zero values
  cens <- as.integer(!is.na(d$cling_leaf) & d$cling_leaf == 0)
  leaf <- dplyr::if_else(cens == 1, NA_real_, log(d$cling_leaf))
  leaf_mean <- tapply(leaf, s, mean, na.rm = TRUE)
  leaf_sd <- tapply(leaf, s, sd, na.rm = TRUE)
  c_leaf <- as.vector((log(leaf_limit) - leaf_mean) / leaf_sd)

  u <- cbind(
    sprint = std_by(log(d$sprint), s),
    bite = std_by(log(d$bite), s),
    bark = std_by(log(d$cling_bark), s),
    leaf = (leaf - leaf_mean[s]) / leaf_sd[s]
  )
  # observed = measured (censored leaf values count as observed)
  obs <- 1L * !is.na(u)
  obs[, "leaf"] <- as.integer(!is.na(d$cling_leaf))
  u[is.na(u)] <- 0 # placeholder; missing/censored values are parameters

  list(
    N = nrow(d),
    S = length(species_levels),
    J = length(perf_traits),
    K = ncol(x),
    species = s,
    male = as.integer(d$sex == "male"),
    x = x,
    w = w,
    shape = match(shape_of, morph_traits),
    u = u,
    obs = obs,
    cens = cens,
    c_leaf = c_leaf,
    y2 = d$seen_fall,
    y3 = d$surv_12m
  )
}
