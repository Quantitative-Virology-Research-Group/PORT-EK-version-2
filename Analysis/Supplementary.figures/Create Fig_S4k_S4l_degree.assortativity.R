Compute.network.assortativity <- function(df) {
  
  df <- df %>% dplyr::filter(sum != 0)
  
  df.btp <- sample_n(df, 1000, replace = T)
  
  df.input <- df %>% dplyr::select(isolate.x, isolate.y, distance)
  
  network.igraph <- graph_from_data_frame(d = df.input, directed = F)
  
  degree.assortativity <- assortativity.degree(network.igraph)
  
  return(degree.assortativity)
}
