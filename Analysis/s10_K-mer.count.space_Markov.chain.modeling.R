## ============================================================================
## REVISED ANALYSIS (answers Reviewer 1, Major 7)
##  * Network: nodes = isolates (n_per_subtype per group); features = marker k-mers; each isolate is joined to its
##    k nearest neighbours (undirected); edge weight = mean over markers of (x_i + x_j)/2.  Zero-weight edges are removed.
##  * Discrete-time Markov chain on the isolates: P_ij = w_ij / deg(i)   (isolates without any edge have an all-zero row).
##  * One-step probabilities p_within / p_across per isolate.
##  * A LABEL-PERMUTATION NULL on the same graph is included.
##  * Long-range structure: expected fraction of the first `walk_steps` steps spent in the starting group, computed EXACTLY from
##    powers of P (no random simulation; an optional Monte-Carlo check is available with n_sim_walks > 0), with a label-permutation null.
##  * Connected components and isolated nodes are reported (the graph need not be connected).
##  * Lumped 7x7 group-level one-step matrix = average, over the isolates of group s, of the probability mass that lands in group t.
##    It is a DESCRIPTIVE one-step summary; its powers are compared with the exact multi-step group probabilities.
##  * Markers were selected with the group labels, so strong within-group connectivity is expected: descriptive, not independent validation.
## ============================================================================
required_pkgs <- c("readr", "dplyr", "tidyr", "stringr", "igraph", "scales", "ggplot2")
to_install <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
library(readr); library(dplyr); library(tidyr); library(stringr); library(igraph); library(scales); library(ggplot2)
```

## 1. Parameters (edit here)

## ---- input files written by `portek find_enriched` (same project, corrected pipeline) -------------------------
input <- list(
  DNA = list(counts = "/path_to_file/15mer_counts_for_classifier.csv",
             stats  = "/path_to_file/enriched_15mers_stats.csv"),
  RNA = list(counts = "/path_to_file/15mer_counts_for_classifier.csv",
             stats  = "/path_to_file/enriched_15mers_stats.csv")
)
panels  <- list(DNA = c(network = "A", bars = "C", heatmap = "E", chain = "G"),   # figure panel letters
                RNA = c(network = "B", bars = "D", heatmap = "F", chain = "H"))
out_dir <- "/path_to_file/Figure4"                      # SVG/CSV output folder
k_label <- 15                             # k used in the panel titles

seed               <- 42    # sampling of isolates and all permutations (report this number in the Methods)
n_per_subtype      <- 50    # isolates sampled per group
n_markers_per_sub  <- 50    # top marker k-mers per group (exclusive k-mers first)
k_neighbors        <- 4     # nearest neighbours per node
walk_steps         <- 50    # steps for the occupancy statistic
n_perm_within      <- 1000  # label permutations for p_within
n_perm_occupancy   <- 200   # label permutations for the 50-step occupancy
n_sim_walks        <- 0     # >0: additional Monte-Carlo check with this many independent walks (e.g. 10000); 0 = skip
exclude_isolated_in_lumping <- FALSE  # FALSE reproduces mean(diag) == mean(p_within); TRUE: average only over isolates with edges (rows sum to 1)
min_edge_prob      <- 0.01  # chain diagram (4G/4H): hide between-group arrows below this probability

## group names: R is the name of the rare-subtype pool in the PORT-EK project; it is shown as "Rare"
alias <- c(R = "Rare")
subtype_colors <- c("A" = "purple3", "B" = "orange2", "C" = "blue", "D" = "green4",
                    "Rare" = "pink2", "CRF01.AE" = "red", "CRF02.AG" = "chocolate4")
subtype_labels <- c("A" = "Subtype A", "B" = "Subtype B", "C" = "Subtype C", "D" = "Subtype D",
                    "Rare" = "Rare subtype", "CRF01.AE" = "CRF01-AE", "CRF02.AG" = "CRF02-AG")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


## 2. Functions

unify <- function(x) ifelse(x %in% names(alias), unname(alias[x]), x)

read_flex <- function(path) {                      # comma- or whitespace-separated
  d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(d) == 1) d <- read.table(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  d
}

select_markers <- function(st, subtype, n) {       # exclusive k-mers first, then highest average count, then RMSE
  pn  <- unique(c(subtype, names(alias)[alias == subtype]))            # e.g. "Rare" and "R"
  cand <- st %>% filter(group %in% paste0(pn, "_enriched"))
  avg_col <- paste0(pn, "_avg"); avg_col <- avg_col[avg_col %in% names(st)][1]
  if (nrow(cand) == 0) stop("no enriched k-mers for group ", subtype, " in the stats table")
  cand %>% arrange(desc(exclusivity == "exclusive"), desc(.data[[avg_col]]), RMSE) %>% slice_head(n = n) %>% pull(kmer)
}

## ---- network + Markov chain --------------------------------------------------------------------------------------
build_chain <- function(counts_path, stats_path) {
  cl <- read_flex(counts_path); st <- read_flex(stats_path)
  names(cl)[1] <- "sample_name"
  if (!"kmer" %in% names(st)) names(st)[1] <- "kmer"
  cl$subtype <- unify(str_extract(cl$sample_name, "^[^_]+"))
  stopifnot(all(cl$subtype %in% names(subtype_colors)))
  subtypes <- names(subtype_colors)[names(subtype_colors) %in% cl$subtype]
  markers <- unique(unlist(lapply(subtypes, function(s) select_markers(st, s, n_markers_per_sub))))
  markers <- markers[markers %in% colnames(cl)]
  cat("markers:", length(markers), "\n")

  set.seed(seed)
  sampled <- cl %>% group_by(subtype) %>% group_modify(~ .x[sample(nrow(.x), min(n_per_subtype, nrow(.x))), ]) %>% ungroup()
  sampled$node_id <- make.unique(str_trunc(sampled$sample_name, 28))
  X <- as.matrix(sampled[, markers]); rownames(X) <- sampled$node_id
  n <- nrow(X); D <- as.matrix(dist(X))

  edges <- vector("list", n * k_neighbors); e <- 0
  for (i in seq_len(n)) {
    others <- setdiff(seq_len(n), i)
    for (j in others[order(D[i, others])[seq_len(min(k_neighbors, n - 1))]]) {
      a <- min(i, j); b <- max(i, j); e <- e + 1
      edges[[e]] <- data.frame(from = rownames(X)[a], to = rownames(X)[b], weight = mean((X[a, ] + X[b, ]) / 2))
    }
  }
  edges <- bind_rows(edges) %>% distinct(from, to, .keep_all = TRUE) %>% filter(weight > 0)   # zero weight = no connection
  nodes <- data.frame(name = sampled$node_id, subtype = sampled$subtype, stringsAsFactors = FALSE)
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  V(g)$color <- subtype_colors[V(g)$subtype]; V(g)$size <- 6
  E(g)$width <- if (ecount(g) > 0) rescale(E(g)$weight, to = c(0.5, 5)) else numeric(0)

  adj <- as.matrix(as_adjacency_matrix(g, attr = "weight", sparse = FALSE))
  deg <- rowSums(adj); iso <- deg == 0
  P <- adj / ifelse(iso, 1, deg)                    # row-stochastic; isolated nodes: all-zero row
  list(g = g, nodes = nodes, P = P, iso = iso, deg = deg, subtypes = subtypes, markers = markers)
}

## ---- statistics --------------------------------------------------------------------------------------------------
analyse_chain <- function(ch) {
  P <- ch$P; sv <- ch$nodes$subtype; n <- length(sv); iso <- ch$iso
  same <- outer(sv, sv, "=="); diag(same) <- FALSE
  p_within <- rowSums(P * same); p_across <- rowSums(P * !same)

  ## label-permutation null for the one-step probability (same graph, labels shuffled over nodes)
  set.seed(seed)
  null_w <- replicate(n_perm_within, { pv <- sample(sv); mean(rowSums(P * { s <- outer(pv, pv, "=="); diag(s) <- FALSE; s })) })
  obs_w  <- mean(p_within)
  p_emp_w <- (1 + sum(null_w >= obs_w)) / (1 + n_perm_within)

  ## exact expected occupancy of the starting group over `walk_steps` steps, uniform start
  S <- matrix(0, n, n); Pt <- diag(n)
  for (t in seq_len(walk_steps)) { Pt <- Pt %*% P; S <- S + Pt }
  S <- S / walk_steps
  home <- outer(sv, sv, "==")
  occ_node <- rowSums(S * home); obs_o <- mean(occ_node)
  null_o <- replicate(n_perm_occupancy, { pv <- sample(sv); mean(rowSums(S * outer(pv, pv, "=="))) })
  p_emp_o <- (1 + sum(null_o >= obs_o)) / (1 + n_perm_occupancy)

  ## connected components of the weighted graph
  comp <- components(ch$g); csize <- comp$csize
  comp_tab <- as.data.frame.matrix(table(component = comp$membership, group = sv))
  comp_tab$size <- csize[as.integer(rownames(comp_tab))]

  ## lumped group-level one-step matrix
  groups <- ch$subtypes
  keep <- if (exclude_isolated_in_lumping) !iso else rep(TRUE, n)
  lump <- function(M) { out <- matrix(0, length(groups), length(groups), dimnames = list(groups, groups))
    for (s in groups) for (t in groups) out[s, t] <- mean(rowSums(M[which(sv == s & keep), which(sv == t), drop = FALSE])); out }
  Pl <- lump(P)
  matpow <- function(M, k) { R <- diag(nrow(M)); for (i in seq_len(k)) R <- R %*% M; R }
  lump_err <- sapply(c(2, 5, 10), function(t) max(abs(matpow(Pl, t) - lump(matpow(P, t)))))
  names(lump_err) <- paste0("t", c(2, 5, 10))

  ## optional Monte-Carlo check of the exact occupancy
  sim_occ <- NA
  if (n_sim_walks > 0) {
    set.seed(seed + 1); cumP <- t(apply(P, 1, cumsum)); has <- rowSums(P) > 0
    cur <- sample(n, n_sim_walks, replace = TRUE); start <- cur; hit <- numeric(n_sim_walks)
    for (s in seq_len(walk_steps)) {
      u <- runif(n_sim_walks); nxt <- max.col(cumP[cur, , drop = FALSE] > u, ties.method = "first")
      cur <- ifelse(has[cur], nxt, cur)               # a walker on an isolated node stays put
      hit <- hit + (sv[cur] == sv[start]) 
    }
    sim_occ <- mean(hit / walk_steps)
  }

  summary <- data.frame(
    nodes = n, markers = length(ch$markers), isolated_nodes = sum(iso), components = length(csize),
    component_sizes = paste(sort(csize, decreasing = TRUE), collapse = ", "),
    mean_p_within = obs_w, mean_p_across = mean(p_across), mean_p_within_excl_isolated = mean(p_within[!iso]),
    perm_null_p_within_mean = mean(null_w), perm_null_p_within_sd = sd(null_w), p_within_empirical_p = p_emp_w,
    occupancy_obs = obs_o, occupancy_obs_excl_isolated = mean(occ_node[!iso]),
    occupancy_null_mean = mean(null_o), occupancy_null_sd = sd(null_o), occupancy_empirical_p = p_emp_o,
    occupancy_monte_carlo = sim_occ, lumped_mean_diag = mean(diag(Pl)),
    lump_err_t2 = lump_err[["t2"]], lump_err_t5 = lump_err[["t5"]], lump_err_t10 = lump_err[["t10"]])
  list(summary = summary, p_within = p_within, p_across = p_across, Pl = Pl, comp_tab = comp_tab,
       stationary_mass = tapply(ch$deg / sum(ch$deg), sv, sum))
}

## ---- figure panels -----------------------------------------------------------------------------------------------
plot_network <- function(ch, file, title, add_labels = TRUE) {            # panels 4A / 4B
  g <- ch$g; set.seed(1)
  lay <- norm_coords(layout_with_fr(g), -1, 1, -1, 1)
  svg(file, width = 7, height = 7); par(mar = c(1, 1, 2, 1))
  plot(g, layout = lay, rescale = FALSE, vertex.label = NA, vertex.color = V(g)$color, vertex.frame.color = "white",
       vertex.size = V(g)$size, edge.width = E(g)$width, edge.color = adjustcolor("gray40", alpha.f = 0.6), main = title)
  if (add_labels) for (s in ch$subtypes) {                       # group name at the median node position, with a white halo
    idx <- which(V(g)$subtype == s); x <- median(lay[idx, 1]); y <- median(lay[idx, 2]); lab <- subtype_labels[[s]]
    for (d in list(c(1, 0), c(-1, 0), c(0, 1), c(0, -1), c(1, 1), c(-1, -1), c(1, -1), c(-1, 1)))
      text(x + d[1] * 0.006, y + d[2] * 0.006, lab, col = "white", cex = 1, font = 2)
    text(x, y, lab, col = subtype_colors[[s]], cex = 1, font = 2) }
  dev.off()
}

plot_bars <- function(ch, an, file, title) {                               # panels 4C / 4D
  df <- data.frame(subtype = ch$nodes$subtype, within = an$p_within, across = an$p_across) %>%
    pivot_longer(c(within, across), names_to = "cond", values_to = "prob") %>%
    group_by(subtype, cond) %>% summarise(prob = mean(prob), .groups = "drop") %>%
    mutate(cond = factor(cond, levels = c("within", "across")), subtype = factor(subtype, levels = ch$subtypes))
  p <- ggplot(df, aes(x = cond, y = prob, fill = cond)) + geom_bar(stat = "identity", color = "black") +
    facet_grid(. ~ subtype) + theme_bw() + scale_fill_manual(values = c(within = "orange", across = "blue")) +
    labs(title = title, x = NULL, y = "One-step transition probability", fill = NULL) +
    theme(axis.text.x = element_text(size = 11, colour = "black", angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 11, colour = "black"), plot.title = element_text(size = 11))
  svg(file, width = 8, height = 4); print(p); dev.off()
}

plot_heatmap <- function(Pl, file, title) {                                # panels 4E / 4F
  ord <- names(subtype_colors)[names(subtype_colors) %in% rownames(Pl)]; M <- Pl[ord, ord]
  svg(file, width = 8, height = 8); par(mar = c(7, 7, 4, 2))
  image(1:ncol(M), 1:nrow(M), t(M)[, nrow(M):1], axes = FALSE, xlab = "", ylab = "", zlim = c(0, 1),
        col = hcl.colors(50, "YlOrRd", rev = TRUE), main = title)
  axis(1, at = 1:ncol(M), labels = subtype_labels[colnames(M)], las = 2, cex.axis = 0.85)
  axis(2, at = 1:nrow(M), labels = rev(subtype_labels[rownames(M)]), las = 2, cex.axis = 0.85)
  mtext("Group at the next step", side = 1, line = 5.5); mtext("Group at the current step", side = 2, line = 5.5)
  for (i in 1:nrow(M)) for (j in 1:ncol(M)) text(j, nrow(M) - i + 1, sprintf("%.2f", M[i, j]), cex = 0.75, col = if (M[i, j] > 0.5) "white" else "black")
  box(); dev.off()
}

plot_chain <- function(Pl, file, title) {                                  # panels 4G / 4H
  ord <- names(subtype_colors)[names(subtype_colors) %in% rownames(Pl)]; M <- Pl[ord, ord]
  M[row(M) != col(M) & M < min_edge_prob] <- 0                             # hide negligible between-group arrows
  gs <- graph_from_adjacency_matrix(M, mode = "directed", weighted = TRUE, diag = TRUE)
  V(gs)$color <- subtype_colors[V(gs)$name]; E(gs)$width <- rescale(E(gs)$weight, to = c(0.5, 10))
  lay <- layout_in_circle(gs); ang <- atan2(lay[, 2], lay[, 1])                # self-loops point outwards
  ee <- ends(gs, E(gs), names = FALSE); loop_ang <- ifelse(ee[, 1] == ee[, 2], ang[ee[, 1]], 0)
  svg(file, width = 8, height = 8)
  plot(gs, layout = lay, xlim = c(-1.5, 1.5), ylim = c(-1.4, 1.4), vertex.label = subtype_labels[V(gs)$name], vertex.label.cex = 0.8, vertex.label.dist = 3,
       vertex.size = 22, edge.width = E(gs)$width, edge.arrow.size = 0.4, edge.curved = 0.25, edge.loop.angle = loop_ang,
       main = title)
  dev.off()
}

## 3. Run DNA and RNA

results <- list()
for (mol in names(input)) {
  cat("\n=====", mol, "=====\n")
  ch <- build_chain(input[[mol]]$counts, input[[mol]]$stats)
  an <- analyse_chain(ch); results[[mol]] <- an$summary

  cat("connected components (rows) x group (columns):\n"); print(an$comp_tab)
  cat("stationary (degree-proportional) mass per group (meaningful only within components):\n"); print(round(an$stationary_mass, 3))
  cat("lumped group-level one-step matrix:\n"); print(round(an$Pl, 3))
  print(t(an$summary))

  pl <- panels[[mol]]; tt <- paste0("Enriched ", mol, " ", k_label, "-mers")
  plot_network(ch, file.path(out_dir, sprintf("Fig4%s_network_%s.svg", pl[["network"]], mol)), tt)
  plot_bars(ch, an, file.path(out_dir, sprintf("Fig4%s_p_within_p_across_%s.svg", pl[["bars"]], mol)), tt)
  plot_heatmap(an$Pl, file.path(out_dir, sprintf("Fig4%s_lumped_matrix_%s.svg", pl[["heatmap"]], mol)), tt)
  plot_chain(an$Pl, file.path(out_dir, sprintf("Fig4%s_lumped_chain_%s.svg", pl[["chain"]], mol)), tt)

  write.csv(an$Pl, file.path(out_dir, sprintf("lumped_matrix_%s.csv", mol)))
  write.csv(data.frame(node = ch$nodes$name, group = ch$nodes$subtype, p_within = an$p_within, p_across = an$p_across, isolated = ch$iso),
            file.path(out_dir, sprintf("node_summary_%s.csv", mol)), row.names = FALSE)
  write.csv(an$comp_tab, file.path(out_dir, sprintf("components_%s.csv", mol)))
}
summary_all <- bind_rows(results, .id = "molecule")
write.csv(summary_all, file.path(out_dir, "random_walk_summary.csv"), row.names = FALSE)
print(t(summary_all))
