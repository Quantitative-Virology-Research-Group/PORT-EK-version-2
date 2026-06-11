kmer.count.matrix <- function(df) {
  df$id <- c(1:nrow(df)) # indexing K-mers
  df.samp <- sample_n(df, 100, replace = T)
  df1 <- df.samp %>% dplyr::select(isolate, kmer.count, subtype, id) %>% dplyr::mutate(temp = "1")
  df2 <- df.samp %>% dplyr::select(isolate, kmer.count, subtype, id) %>% dplyr::mutate(temp = "1")
  
  df.merge <- dplyr::full_join(df1, df2, by = "temp", relationship = "many-to-many")
  
df.fi <- data.frame()
  
  for(i in 1:nrow(df.merge)) {
      if (isTRUE(df.merge[i,1] == df.merge[i,6])) {
        next
      }
      df.merge.tmp <- as.data.frame(df.merge[i,1]) %>% dplyr::mutate(id.x = df.merge[i,4], ,kmer.x = df.merge[i,2], isolate.y = df.merge[i,6], id.y = df.merge[i,9], kmer.y = df.merge[i,7], sum = df.merge[i,2]+df.merge[i,7])
      names(df.merge.tmp) <- c("isolate.x", "id.x","kmer.x", "isolate.y", "id.y", "kmer.y","sum")
      df.fi <- rbind(df.fi, df.merge.tmp)
  }
  return(df.fi)
}

Execution.kmer.count.matrix <- function(df, n.subsets = 1000) {
  
  df.output <- purrr::map_dfr(
    .x = seq_len(n.subsets),
    .f = ~ kmer.count.matrix(df) %>% dplyr::mutate(subset = .x),
    .progress = TRUE
  )
  
  return(df.output)
}

Execution.kmer.count.matrix.into.isolate <- function(df) {

  # Build raw pairwise table
  df.raw <- Execution.kmer.count.matrix(df) %>%
    dplyr::mutate(distance = 1- sum)# original distance kept
  
  df.raw <- df.raw %>% dplyr::arrange(distance)
  
  n.top10 <- 0.1*nrow(df.raw)
  n.top20 <- 0.2*nrow(df.raw)
  n.top30 <- 0.3*nrow(df.raw)
  
  df.raw.top10 <- df.raw[1:n.top10,] %>% dplyr::mutate(rank = "top10")
  df.raw.top20 <- df.raw[1:n.top20,] %>% dplyr::mutate(rank = "top20")
  df.raw.top30 <- df.raw[1:n.top30,] %>% dplyr::mutate(rank = "top30")
  
  df.output <- dplyr::bind_rows(df.raw.top10, df.raw.top20, df.raw.top30)
  
  return(df.output)
}

plot.network.isolate <- function(df.pool) {
  
  df.input <- sample_n(df, 1000) %>% dplyr::select(isolate.x, isolate.y, distance)
  
  df.isolate.x <- df %>% dplyr::select(isolate.x) %>% dplyr::rename(isolate = isolate.x)
  df.isolate.y <- df %>% dplyr::select(isolate.y) %>% dplyr::rename(isolate = isolate.y)
  node.isolate <- dplyr::bind_rows(df.isolate.x, df.isolate.y) %>% unique() 
  
  network.igraph <- graph_from_data_frame(d = df.input, vertices = node.isolate, directed = F)
  
  plot(network.igraph, vertex.size = 10, vertex.label = NA)
}
