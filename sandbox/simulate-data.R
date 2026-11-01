
## Packages ----

library(tidyverse)

# Scratch space. Shared code lives in R/; see 03-validation.qmd for the
# prior predictive checks and simulation-based calibration.

source("R/prepare-data.R")
source("R/priors.R")
source("R/simulate.R")

## Test ----

lizards <- read_csv("data/lizards-spring-2024.csv", show_col_types = FALSE)
design <- prepare_data(lizards)

set.seed(1)
draw <- simulate_one(design, priors, direct_effect = TRUE)
length(draw$variables)
