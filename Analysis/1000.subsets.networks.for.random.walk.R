##====================================================================================##
# functions
# 1. Generation of 1000 graph matrices (100 nodes in each network)
isolate.graph.matrix <- function(df) {
  df$id <- c(1:nrow(df)) # indexing isolates
  df.samp <- sample_n(df, 100, replace = T)
  df1 <- df.samp %>% dplyr::select(isolate, kmer.count, subtype, id) %>% dplyr::mutate(temp = "1")
  df2 <- df.samp %>% dplyr::select(isolate, kmer.count, subtype, id) %>% dplyr::mutate(temp = "1")
  
  df.merge <- dplyr::full_join(df1, df2, by = "temp", relationship = "many-to-many")
  
df.fi <- data.frame()
  
  for(i in 1:nrow(df.merge)) {
      if (isTRUE(df.merge[i,1] == df.merge[i,6])) {
        next
      }
      df.merge.tmp <- as.data.frame(df.merge[i,1]) %>% dplyr::mutate(id.x = df.merge[i,4], isolate.y = df.merge[i,6], id.y = df.merge[i,9], sum = df.merge[i,2] + df.merge[i,7])
      names(df.merge.tmp) <- c("isolate.x", "id.x", "isolate.y", "id.y","sum")
      df.fi <- rbind(df.fi, df.merge.tmp)
  }
  return(df.fi)
}

Execution.kmer.count.matrix.into.isolate <- function(df) {
  df.A <- df %>% dplyr::filter(subtype == "A")
  df.B <- df %>% dplyr::filter(subtype == "B")
  df.C <- df %>% dplyr::filter(subtype == "C")
  df.D <- df %>% dplyr::filter(subtype == "D")
  df.rest <- df %>% dplyr::filter(subtype == "rest")
  
  df.A.output <- isolate.graph.matrix(df.A) %>% dplyr::mutate(subtype = "A")
  df.B.output <- isolate.graph.matrix(df.B) %>% dplyr::mutate(subtype = "B")
  df.C.output <- isolate.graph.matrix(df.C) %>% dplyr::mutate(subtype = "C")
  df.D.output <- isolate.graph.matrix(df.D) %>% dplyr::mutate(subtype = "D")
  df.rest.output <- isolate.graph.matrix(df.rest) %>% dplyr::mutate(subtype = "rest")

  df.isolate.output <- dplyr::bind_rows(df.A.output, df.B.output, df.C.output, df.D.output, df.rest.output)
  
  return(df.isolate.output)
}

# 2. reshape graph matrix
reshape.graph.matrix <- function(df) {
  df_ <- df %>% dplyr::select(id.x, id.y, sum) %>% filter(sum != 0)
  df.complete <- df_ %>% complete(id.x, id.y, fill = list(sum = 0))
  df.wide <- dcast(df.complete, id.x ~ id.y, value.var = "sum", fun.aggregate = sum)
  rownames(df.wide) <- df.wide$id.x
  df.wide$id.x <- NULL
  output.mx <- as.matrix(df.wide)
  
  return(output.mx)
}

## 3. A random walk of 10k steps and calculate the probability
### ── a. RANDOM WALK SIMULATOR ──────────────────────────────────────────────────
simulate_random_walk.subtype <- function(graph_matrix, start_node, num_steps) {
  path <- integer(num_steps)
  path[1]     <- start_node
  current_node <- start_node

  for (i in 2:num_steps) {
    neighbors <- which(graph_matrix[current_node, ] != 0)
    if (length(neighbors) == 0) {          # dead-end: stop early
      path <- path[1:(i - 1)]
      break
    }
    current_node <- sample(neighbors, 1)   # unweighted step
    path[i]      <- current_node
  }
  return(path)
}

### ── b. HELPERS: node-name → subtype / start-node detection ───────────────────
#' Map a single node name to its subtype label
node_to_subtype <- function(name) {
  dplyr::case_when(
    startsWith(name, "A_")    ~ "A",
    startsWith(name, "B_")    ~ "B",
    startsWith(name, "C_")    ~ "C",
    startsWith(name, "D_")    ~ "D",
    startsWith(name, "rest_") ~ "R",
    TRUE                      ~ NA_character_
  )
}

#' Return the index of the first node whose name matches a known prefix.
#' Falls back to node 1 if none found.
get_start_node <- function(graph_matrix) {
  nms <- rownames(graph_matrix)
  if (is.null(nms)) return(1L)

  idx <- which(grepl("^(A_|B_|C_|D_|rest_)", nms))
  if (length(idx) > 0L) idx[1L] else 1L
}

### ── c. PROBABILITY CALCULATOR ─────────────────────────────────────────────────
#' Compute visit-probability per subtype and same-scenario return probability.
#'
#' @param path        Integer vector of visited node indices.
#' @param node_names  Character vector – rownames of the graph matrix.
#' @return data.frame with columns: subtype, prob, start_subtype, same_scenario_prob
random.walk.prob <- function(path, node_names) {
  visited_subtypes <- node_to_subtype(node_names[path])
  #start_subtype    <- visited_subtypes[1]
  start_subtype <- visited_subtypes[[sample(1:length(visited_subtypes), 1)]]
  total            <- length(path)

  all_subtypes <- c("A", "B", "C", "D", "R")

  prob_df <- data.frame(
    subtype            = all_subtypes,
    prob               = sapply(all_subtypes,
                                function(s) sum(visited_subtypes == s,
                                                na.rm = TRUE) / total),
    start_subtype      = start_subtype,
    same_scenario_prob = sum(visited_subtypes == start_subtype,
                             na.rm = TRUE) / total,
    stringsAsFactors   = FALSE
  )
  return(prob_df)
}

### ── 4. MAIN WRAPPER – works for any list length ───────────────────────────────
#' Run random walks on every graph in graph.list and return probability tables.
#'
#' @param graph.list  List of adjacency matrices (any length, e.g. 1000).
#'                    Each matrix MUST have rownames following the
#'                    "A_", "B_", "C_", "D_", "rest_" naming convention.
#' @param num_steps   Number of walk steps (default 10 000).
#' @return List of data.frames (one per graph), each with walk statistics.
walk.path.list <- function(graph.list, num_steps = 10000L) {

  results <- lapply(seq_along(graph.list), function(i) {

    g          <- graph.list[[i]]
    node_names <- rownames(g)

    # Auto-detect start node from node names
    start_node  <- get_start_node(g)
    start_label <- if (!is.null(node_names)) node_names[start_node] else NA

    # Simulate
    path <- simulate_random_walk.subtype(g, start_node, num_steps)

    # Compute probabilities
    prob_df <- random.walk.prob(path, node_names)

    # Attach metadata
    prob_df$graph_id         <- i
    prob_df$start_node_index <- start_node
    prob_df$start_node_name  <- start_label
    prob_df$actual_steps     <- length(path)   # may be < num_steps if dead-end

    return(prob_df)
  })

  # Optionally combine into a single flat data.frame for easy analysis
  all_results <- dplyr::bind_rows(results)
  return(all_results)
}
##====================================================================================##

#execution
isolate.kmer.count.mx.list <- parallel::mclapply(1:1000, function(i) {
    Execution.kmer.count.matrix.into.isolate(isolate.kmer.count) %>% dplyr::mutate(subset = as.character(i))
  }, mc.cores = parallel::detectCores() - 1)

graph.matrix.isolate.list <- lapply(1:1000, function(i) {
  reshape.graph.matrix.isolate(isolate.kmer.count.mx.list[[i]])
})

df.prob <- walk.path.list(graph.matrix.isolate.list, num_steps = 10000)

