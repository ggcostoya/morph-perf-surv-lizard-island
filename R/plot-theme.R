# Shared plot styling. Species colours are the first three slots of a
# colour-blind-safe categorical palette (validated for all pairs).

species_colours <- c(caro = "#2a78d6", cris = "#eb6834", sag = "#1baf7a")
species_names <- c(
  caro = "A. carolinensis",
  cris = "A. cristatellus",
  sag = "A. sagrei"
)

scale_colour_species <- function(...) {
  ggplot2::scale_colour_manual(values = species_colours, labels = species_names, name = NULL, ...)
}
scale_fill_species <- function(...) {
  ggplot2::scale_fill_manual(values = species_colours, labels = species_names, name = NULL, ...)
}

theme_lizard <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#ecebe8", linewidth = 0.3),
      legend.position = "top",
      legend.justification = "left",
      legend.text = ggplot2::element_text(face = "italic"),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0),
      axis.title = ggplot2::element_text(colour = "#52514e")
    )
}
