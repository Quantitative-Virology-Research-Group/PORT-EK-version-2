## ============================================================================
## HIV-1 Subtype Spatial-Distribution Network + Markov Chain Random-Walk Model
## ----------------------------------------------------------------------------
## PART 1 (Sections 1-6)  builds the same weighted, undirected isolate network
##                         as before: nodes = isolates, colored by subtype;
##                         edges = k-NN in marker-k-mer space, weight = average
##                         k-mer count between the two connected isolates.
##
## PART 2 (Sections 7-10) treats a random walk moving over that weighted
##                         network as a discrete-time MARKOV chain.
##                        This is done by:
##                           7.  building the random-walk transition matrix P,
##                           8.  computing, per node, the one-step probability
##                               of remaining within vs. leaving its subtype,
##                               and testing the difference,
##                          9.  lumping the isolate-level chain into a 7x7
##                               subtype-level Markov chain,
##                          10.  validating with a Monte-Carlo simulation of
##                               multi-step random walks.
## ============================================================================
 
## ---- 0. Packages -----------------------------------------------------------
required_pkgs <- c("readr", "dplyr", "tidyr", "stringr", "igraph", "scales")
to_install <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
 
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(igraph)
library(scales)
 
set.seed(42)
 
## ---- 1. User-tunable parameters --------------------------------------------
n_per_subtype      <- 50   # isolates sampled per subtype (keeps the graph readable)
n_markers_per_sub  <- 50   # top marker k-mers used per subtype to build feature space
k_neighbors        <- 4    # number of nearest neighbours connected per node
 
walk_steps         <- 50   # steps per simulated random walk (Section 11)
walk_n_starts      <- 10000  # number of independent walks simulated (Section 11)
 
## Official subtype color code
subtype_colors <- c(
  "A"          = "purple3",
  "B"          = "orange2",
  "C"          = "blue",
  "D"          = "green4",
  "Rare"       = "pink2",
  "CRF01.AE"   = "red",
  "CRF02.AG"   = "chocolate4"
)
subtype_labels <- c(
  "A"        = "Subtype A",
  "B"        = "Subtype B",
  "C"        = "Subtype C",
  "D"        = "Subtype D",
  "Rare"     = "Rare subtype",
  "CRF01.AE" = "CRF01-AE",
  "CRF02.AG" = "CRF02-AG"
)
 
## ---- 2. Load data -----------------------------------------------------------
classifier <- read.table("path_to_file", header = T, stringsAsFactors = F) # isolate x k-mer matrix
kmer_info  <- read.table("path_to_file", header = T, stringsAsFactors = F) # k-mer -> subtype enrichment metadata
 
classifier <- classifier %>%
  mutate(subtype = str_extract(sample_name, "^[^_]+")) %>%
  relocate(subtype, .after = sample_name)
 
stopifnot(all(classifier$subtype %in% names(subtype_colors)))
 
## ---- 3. Select subtype-specific marker k-mers -------------------------------
select_markers <- function(subtype, n) {
  grp <- paste0(subtype, "_enriched")
  avg_col <- paste0(subtype, "_avg")
 
  candidates <- kmer_info %>% filter(group == grp)
 
  candidates %>%
    arrange(desc(exclusivity == "exclusive"), desc(.data[[avg_col]]), RMSE) %>%
    slice_head(n = n) %>%
    pull(kmer)
}
 
subtypes <- names(subtype_colors)
marker_list <- setNames(
  lapply(subtypes, select_markers, n = n_markers_per_sub),
  subtypes
)
marker_kmers <- unique(unlist(marker_list))
marker_kmers <- marker_kmers[marker_kmers %in% colnames(classifier)]
 
## ---- 4. Sample isolates per subtype ------------
sampled <- classifier %>%
  group_by(subtype) %>%
  slice_sample(n = n_per_subtype) %>%
  ungroup()
 
sampled <- sampled %>%
  mutate(node_id = make.unique(str_trunc(sample_name, 28)))
 
feature_mat <- as.matrix(sampled[, marker_kmers])
rownames(feature_mat) <- sampled$node_id
 
## ---- 5. Build a weighted, undirected k-nearest-neighbour graph -------------
dist_mat <- as.matrix(dist(feature_mat, method = "euclidean"))
n_nodes  <- nrow(feature_mat)
edge_rows <- list()
 
for (i in seq_len(n_nodes)) {
  nn_idx <- order(dist_mat[i, -i])[seq_len(min(k_neighbors, n_nodes - 1))]
  nn_idx <- setdiff(seq_len(n_nodes), i)[nn_idx]
  for (j in nn_idx) {
    a <- min(i, j); b <- max(i, j)
    ## average k-mer count between the two adjacent nodes across the marker
    ## k-mer feature space (mean of the per-kmer average of the two isolates)
    weight <- mean((feature_mat[a, ] + feature_mat[b, ]) / 2)
    edge_rows[[length(edge_rows) + 1]] <- data.frame(
      from = rownames(feature_mat)[a],
      to   = rownames(feature_mat)[b],
      weight = weight
    )
  }
}
 
edges_df <- bind_rows(edge_rows) %>% distinct(from, to, .keep_all = TRUE)
 
## ---- 6. Assemble the igraph object & plot the base network -----------------
nodes_df <- sampled %>% transmute(name = node_id, subtype = subtype)
 
g <- graph_from_data_frame(edges_df, directed = FALSE, vertices = nodes_df)
V(g)$color <- subtype_colors[V(g)$subtype]
V(g)$size  <- 6
E(g)$width <- rescale(E(g)$weight, to = c(0.5, 5))
 
set.seed(1)
layout_xy <- layout_with_fr(g)
 
plot(
  g,
  layout             = layout_xy,
  vertex.label       = NA,
  vertex.color       = V(g)$color,
  vertex.frame.color = "white",
  vertex.size        = V(g)$size,
  edge.width         = E(g)$width,
  edge.color         = adjustcolor("gray40", alpha.f = 0.6),
  main               = "Spatial Distribution Network of HIV-1 Isolates Across Subtypes"
)
legend(
  "topright",
  legend = subtype_labels[names(subtype_colors)],
  col    = subtype_colors,
  pch    = 19, pt.cex = 1.4, bty = "n", cex = 0.8,
  title  = "HIV-1 Subtype"
)
 
## ============================================================================
## PART 2 -- MARKOV CHAIN / RANDOM-WALK MODEL
## ============================================================================
 
## ---- 7. Random-walk transition matrix (node level) -------------------------
adj <- as.matrix(as_adjacency_matrix(g, attr = "weight", sparse = FALSE))
deg <- rowSums(adj)
deg[deg == 0] <- 1e-9          # guard against isolated nodes
P <- adj / deg
 
stopifnot(all(abs(rowSums(P) - 1) < 1e-6))   # P must be row-stochastic
  
## ---- 8. Node-level within- vs. across-subtype one-step probability --------
subtype_vec <- nodes_df$subtype
same_subtype_mat <- outer(subtype_vec, subtype_vec, "==")
diag(same_subtype_mat) <- FALSE   # "moving to itself" isn't a real transition here
 
p_within <- rowSums(P * same_subtype_mat)
p_across <- rowSums(P * (!same_subtype_mat))
## (p_within + p_across == 1 for every node, since every outgoing step is
##  either same-subtype or different-subtype)
 
node_summary <- tibble(
  node             = nodes_df$name,
  subtype          = nodes_df$subtype,
  p_within_subtype = p_within,
  p_across_subtype = p_across
)
 
cat("\nMean one-step probability of staying WITHIN the starting subtype:",
    round(mean(p_within), 3), "\n")
cat("Mean one-step probability of jumping ACROSS to another subtype:   ",
    round(mean(p_across), 3), "\n")
cat("Ratio (within / across):", round(mean(p_within) / mean(p_across), 2), "\n")
 
wilcox_result <- wilcox.test(node_summary$p_within_subtype,
                              node_summary$p_across_subtype,
                              paired = TRUE)
cat("\nWilcoxon signed-rank test (paired by isolate), H1: within > across:\n")
print(wilcox_result)
 
## ---- 9. Lump the isolate-level chain into a 7x7 subtype-level chain ------
## For subtypes s and t, the lumped transition probability is the average,
## over isolates currently in s, of their total outgoing probability mass
## landing on any isolate in t. This yields a smaller, interpretable
## subtype x subtype Markov chain.
subtypes_present <- sort(unique(nodes_df$subtype))
P_subtype <- matrix(0, nrow = length(subtypes_present), ncol = length(subtypes_present),
                     dimnames = list(subtypes_present, subtypes_present))
 
for (s in subtypes_present) {
  idx_s <- which(nodes_df$subtype == s)
  for (t in subtypes_present) {
    idx_t <- which(nodes_df$subtype == t)
    P_subtype[s, t] <- mean(rowSums(P[idx_s, idx_t, drop = FALSE]))
  }
}
P_subtype <- P_subtype / rowSums(P_subtype)   # renormalize row-stochastic
 
cat("\nLumped 7x7 subtype-level random-walk transition matrix:\n")
print(round(P_subtype, 3))
 
within_vals <- diag(P_subtype)
across_vals <- P_subtype[row(P_subtype) != col(P_subtype)]
cat("\nSubtype-level summary:\n")
cat("  Mean diagonal (within-subtype) probability :", round(mean(within_vals), 3), "\n")
cat("  Mean off-diagonal (across-subtype) probability:", round(mean(across_vals), 3), "\n")
cat("  Ratio:", round(mean(within_vals) / mean(across_vals), 2), "\n")
 
## ---- 10. Monte-Carlo validation: simulate multi-step random walks ---------
subtype_of <- setNames(nodes_df$subtype, nodes_df$name)
 
simulate_walk <- function(start, steps, P) {
  current <- start
  home_subtype <- subtype_of[[start]]
  same_steps <- 0L
  for (s in seq_len(steps)) {
    current <- sample(colnames(P), size = 1, prob = P[current, ])
    if (subtype_of[[current]] == home_subtype) same_steps <- same_steps + 1L
  }
  same_steps / steps
}
 
set.seed(123)
starts <- sample(rownames(P), walk_n_starts, replace = TRUE)
sim_fraction_same <- sapply(starts, simulate_walk, steps = walk_steps, P = P)
 
cat("\nMonte-Carlo simulation (", walk_n_starts, " walks x ", walk_steps, " steps):\n", sep = "")
cat("  Mean fraction of time spent within the STARTING subtype:",
    round(mean(sim_fraction_same), 3), "\n")
cat("  Expected fraction under uniform mixing (no subtype structure): ~",
    round(1 / length(subtypes_present), 3), "\n")
