# morph-perf-surv-lizard-island

Does morphology affect survival only through performance? An open,
reproducible analysis of natural selection in three *Anolis* species
on Lizard Island, written as a Quarto website.

- `index.qmd`, `01-data.qmd` to `05-synthesis.qmd`: the pages of the site
- `data/`: cleaned analysis dataset
- `R/`: shared functions (data processing, priors, simulator, plot theme)
- `figures/`: static figures

Render with `quarto render`. The raw field data are not included;
`R/process-raw-data.R` documents how the cleaned dataset was produced
from them.
