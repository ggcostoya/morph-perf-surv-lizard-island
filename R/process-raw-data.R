library(tidyverse)

# Builds the public analysis dataset from the raw field data. The raw file
# (raw_data/, not in the public repository) is read only here. Outputs:
#   data/lizards-spring-2024.csv  one row per spring 2024 lizard of the three
#                                 focal species, cleaned, with flags
#   data/capture-summary.csv      captures per season and species

## Read ----

# "#N/A" (Excel) and "na" (sex, mass) are also missing values
raw <- read.csv(
  "raw_data/combined_lizard_island_GGC.csv",
  na.strings = c("NA", "", "#N/A", "na")
) |>
  dplyr::filter(!dplyr::if_all(dplyr::everything(), is.na)) |>
  dplyr::select(-X) |>
  as_tibble()

## Fall 2024 recapture census ----

# fall 2024 only records captures (no morphology or performance); used
# here to flag which spring lizards were detected at the 6-month census
fall <- raw |>
  dplyr::filter(sampling.season == "fall.2024") |>
  dplyr::mutate(date = as.Date(date, format = "%m/%d/%Y"))
fall_tags <- fall$alpha_tag

# sex at the latest recapture; older/larger lizards are easier to sex
fall_sex <- fall |>
  dplyr::filter(!is.na(sex)) |>
  dplyr::arrange(alpha_tag, date, time) |>
  dplyr::group_by(alpha_tag) |>
  dplyr::summarise(sex_fall = dplyr::last(sex))

## Capture summary (all species and seasons) ----

spring_tags <- raw |>
  dplyr::filter(sampling.season == "spring.2024") |>
  dplyr::pull(alpha_tag)

capture_summary <- raw |>
  dplyr::mutate(
    season = dplyr::recode(sampling.season, spring.2024 = "spring 2024", fall.2024 = "fall 2024"),
    recapture = sampling.season == "fall.2024" & alpha_tag %in% spring_tags
  ) |>
  dplyr::group_by(season, species) |>
  dplyr::summarise(
    captures = dplyr::n(),
    recaptures_of_spring_lizards = sum(recapture),
    .groups = "drop"
  )

## Spring 2024 cohort ----

data <- raw |>
  dplyr::filter(
    sampling.season == "spring.2024",
    species %in% c("carolinensis", "cristatellus", "sagrei")
  ) |>
  dplyr::select(
    lizard_id = lizard.id,
    alpha_tag,
    date,
    species,
    sex,
    mass = mass_g,
    svl,
    head_length = head.length,
    head_width = head.width,
    pectoral_width = pectoral.width,
    upper_forelimb = upper.forelimb,
    lower_forelimb = lower.forelimb,
    forelimb, # = upper_forelimb + lower_forelimb
    forefoot,
    pelvic_width = pelvic.width,
    upper_hindlimb = upper.hindlimb,
    lower_hindlimb = lower.hindlimb,
    hindlimb, # = upper_hindlimb + lower_hindlimb
    hindfoot,
    hind_toe = hind.toe,
    forefoot_width,
    forefoot_area,
    hindfoot_width,
    hindfoot_area,
    sprint = max_sprint_1_mps, # best trial; values are cm/s, not m/s
    tb_sprint = sprint_body_temp,
    bite = max_bite,
    tb_bite = bite_body_temp,
    cling_bark = max_cling_bark,
    cling_leaf = max_cling_leaf,
    surv_6m = survival.6month,
    surv_12m = survival.12month
  ) |>
  dplyr::mutate(
    species = dplyr::recode(
      species,
      carolinensis = "caro", cristatellus = "cris", sagrei = "sag"
    ),
    species = factor(species, levels = c("caro", "cris", "sag")),
    date = as.Date(date, format = "%m/%d/%Y"),
    toepad_area = forefoot_area + hindfoot_area,
    seen_fall = as.integer(alpha_tag %in% fall_tags)
  )

## Corrections ----

# each correction is recorded in `corrections` so it can be traced
is_true <- function(condition) !is.na(condition) & condition
note <- function(notes, condition, text) {
  condition <- is_true(condition)
  dplyr::if_else(condition, paste(notes, text, sep = "; "), notes)
}

data <- data |>
  dplyr::left_join(fall_sex, by = "alpha_tag") |>
  dplyr::group_by(species) |>
  dplyr::mutate(
    hind_toe_z = (log(hind_toe / svl) - median(log(hind_toe / svl))) / mad(log(hind_toe / svl))
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    corrections = "",
    # sex at the latest record (fall recapture) where the two disagree
    corrections = note(corrections, sex_fall != sex, "sex from fall recapture"),
    sex = factor(dplyr::coalesce(sex_fall, sex), levels = c("female", "male")),
    # 4 February-March dates are most likely typos (season ran May-June)
    corrections = note(corrections, format(date, "%m") %in% c("02", "03"), "implausible date removed"),
    date = dplyr::if_else(is_true(format(date, "%m") %in% c("02", "03")), as.Date(NA), date),
    # lizard 714 has its ID entered as sprint body temperature
    corrections = note(corrections, tb_sprint > 45, "sprint temperature was lizard ID, removed"),
    tb_sprint = dplyr::if_else(is_true(tb_sprint > 45), NA_real_, tb_sprint),
    # lizard 402 bite temperature entered as 20.2; trials were run near 30
    corrections = note(corrections, lizard_id == 402 & tb_bite == 20.2, "bite temperature 20.2 corrected to 30.2"),
    tb_bite = dplyr::if_else(is_true(lizard_id == 402 & tb_bite == 20.2), 30.2, tb_bite),
    # lizard 434 scored 0 on both cling surfaces: a failed trial on both
    corrections = note(corrections, cling_bark == 0, "failed cling trial removed"),
    cling_leaf = dplyr::if_else(is_true(cling_bark == 0 & cling_leaf == 0), NA_real_, cling_leaf),
    cling_bark = dplyr::if_else(is_true(cling_bark == 0), NA_real_, cling_bark),
    # remaining leaf zeros are below the 0.01 resolution (left-censored)
    cling_leaf_censored = as.integer(cling_leaf == 0),
    # implausibly short hind toes (robust z of log toe/SVL < -5 within
    # species; clear gap between -7.4 and -4.0) are set to missing
    hind_toe_recorded = hind_toe,
    corrections = note(corrections, hind_toe_z < -5, "implausibly short hind toe removed"),
    hind_toe = dplyr::if_else(is_true(hind_toe_z < -5), NA_real_, hind_toe),
    corrections = dplyr::na_if(sub("^; ", "", corrections), "")
  ) |>
  dplyr::select(-sex_fall, -hind_toe_z)

## Exclusions ----

# lizards are kept in the file but flagged; the analysis uses
# is.na(exclusion). 5 lizards have no ID/tag ("yellowNA") so could never
# be resighted; 2 tags are assigned to two different lizards each;
# morphology is conditioned on in the model, so it must be complete
shared_tag <- duplicated(data$alpha_tag) | duplicated(data$alpha_tag, fromLast = TRUE)
data <- data |>
  dplyr::mutate(
    exclusion = dplyr::case_when(
      is.na(lizard_id) ~ "no tag",
      shared_tag ~ "tag shared by two lizards",
      dplyr::if_any(c(svl, hindlimb, head_width, hind_toe, toepad_area), is.na) ~ "missing morphology"
    ),
    seen_fall = dplyr::if_else(exclusion %in% c("no tag", "tag shared by two lizards"), NA_integer_, seen_fall)
  )

## Checks ----

tracked <- data |> dplyr::filter(!exclusion %in% c("no tag", "tag shared by two lizards"))
stopifnot(
  !anyDuplicated(tracked$alpha_tag),
  all(tracked$surv_12m <= tracked$surv_6m), # alive at 12 months implies alive at 6
  all(tracked$surv_6m >= tracked$seen_fall) # seen in fall implies alive at 6
)

## Save ----

# tag codes are only needed to link recaptures, so they are not released
dir.create("data", showWarnings = FALSE)
data |>
  dplyr::select(-alpha_tag) |>
  dplyr::relocate(exclusion, corrections, .after = sex) |>
  write_csv("data/lizards-spring-2024.csv", na = "NA")
write_csv(capture_summary, "data/capture-summary.csv")
