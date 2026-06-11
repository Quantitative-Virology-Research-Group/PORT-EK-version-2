{r FUN::computation of the "distance-based" similarity matrix}
Compute.distance.based.similarity <- function(df.A, df.B, df.C, df.D, df.R) {
  
  library(ggplot2)
  library(reshape2)
  
  df_ <- dplyr::bind_rows(df.A, df.B, df.C, df.D, df.R) %>% dplyr::select(assort.deg, subtype)

  df <- aggregate(assort.deg ~ subtype, df_, mean)
  
  # Build a pairwise correlation matrix between subtypes
  # We treat each subtype's assort.deg value as a single scalar, so we compute
  # the "distance-based" similarity matrix:  corr(i,j) = 1 - |val_i - val_j| / range
  # Also compute a simple Pearson-style similarity:  corr = 1 - (|diff| / max_diff)
  vals     <- setNames(df$assort.deg, df$subtype)
  subtypes <- df$subtype
  n        <- length(subtypes)
  range_v  <- diff(range(vals))
 
  mat <- outer(vals, vals, FUN = function(x, y) 1 - abs(x - y) / range_v)
  rownames(mat) <- colnames(mat) <- subtypes
  
  # ---- ggplot heatmap ----
  mat_melt <- melt(mat)
  colnames(mat_melt) <- c("Subtype_X", "Subtype_Y", "Similarity")
 
  ggplot(mat_melt, aes(x = Subtype_X, y = Subtype_Y, fill = Similarity)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = round(Similarity, 3)), size = 5, fontface = "bold",
            colour = ifelse(mat_melt$Similarity > 0.85, "white", "black")) +
  scale_fill_gradient2(
    low      = "#4575b4",
    mid      = "#fee090",
    high     = "orangered",
    midpoint = 0.85,
    limits   = c(0.7, 1),
    name     = "Similarity"
  ) +
  labs(
    title    = "Pairwise Similarity Matrix of Subtypes\nbased on Assortativity Degree",
    subtitle = "Similarity = 1 − |assort.deg_i − assort.deg_j| / range",
    x        = "Subtype",
    y        = "Subtype"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 15),
    plot.subtitle   = element_text(hjust = 0.5, colour = "black", size = 11),
    axis.text       = element_text(size = 13, face = "bold"),
    legend.position = "right",
    panel.grid      = element_blank()
  )
}
