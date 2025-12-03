Compute.network.assortativity <- function(df) {
  df <- df %>% dplyr::filter(sum != 0)
  
  df.input <- df %>% dplyr::select(kmer.x, kmer.y, distance)
  
  network.igraph <- graph_from_data_frame(d = df.input, directed = F)
  
  degree.assortativity <- assortativity.degree(network.igraph)
  
  return(degree.assortativity)
}

Compute.network.assortativity.subset <- function(df) {
  df.subset.1 <- df %>% dplyr::filter(subset == "1")
  dg.assort.1 <- Compute.network.assortativity(df.subset.1)
  
  df.subset.2 <- df %>% dplyr::filter(subset == "2")
  dg.assort.2 <- Compute.network.assortativity(df.subset.2)
  
  df.subset.3 <- df %>% dplyr::filter(subset == "3")
  dg.assort.3 <- Compute.network.assortativity(df.subset.3)
  
  df.subset.4 <- df %>% dplyr::filter(subset == "4")
  dg.assort.4 <- Compute.network.assortativity(df.subset.4)
  
  df.subset.5 <- df %>% dplyr::filter(subset == "5")
  dg.assort.5 <- Compute.network.assortativity(df.subset.5)
  
  df.subset.6 <- df %>% dplyr::filter(subset == "6")
  dg.assort.6 <- Compute.network.assortativity(df.subset.6)
  
  df.subset.7 <- df %>% dplyr::filter(subset == "7")
  dg.assort.7 <- Compute.network.assortativity(df.subset.7)
  
  df.subset.8 <- df %>% dplyr::filter(subset == "8")
  dg.assort.8 <- Compute.network.assortativity(df.subset.8)
  
  df.subset.9 <- df %>% dplyr::filter(subset == "9")
  dg.assort.9 <- Compute.network.assortativity(df.subset.9)
  
  df.subset.10 <- df %>% dplyr::filter(subset == "10")
  dg.assort.10 <- Compute.network.assortativity(df.subset.10)
  
  values <- c(dg.assort.1, dg.assort.2, dg.assort.3, dg.assort.4, dg.assort.5, dg.assort.6, dg.assort.7, dg.assort.8, dg.assort.9, dg.assort.10)
  subset <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
  
  df.pool <- as.data.frame(values, subset) %>% na.omit()
  
  return(df.pool)
}

Compute.network.assortativity.subset.subtype <- function(df) {
  
  df.A <- df %>% dplyr::filter(subtype == "A")
  df.A.dg.assort <- Compute.network.assortativity.subset(df.A) %>% dplyr::mutate(subtype = "A")
  
  df.B <- df %>% dplyr::filter(subtype == "B")
  df.B.dg.assort <- Compute.network.assortativity.subset(df.B) %>% dplyr::mutate(subtype = "B")
  
  df.C <- df %>% dplyr::filter(subtype == "C")
  df.C.dg.assort <- Compute.network.assortativity.subset(df.C) %>% dplyr::mutate(subtype = "C")
  
  df.D <- df %>% dplyr::filter(subtype == "D")
  df.D.dg.assort <- Compute.network.assortativity.subset(df.D) %>% dplyr::mutate(subtype = "D")
  
  df.rest <- df %>% dplyr::filter(subtype == "rest")
  df.rest.dg.assort <- Compute.network.assortativity.subset(df.rest) %>% dplyr::mutate(subtype = "rest")
  
  df.fi <- bind_rows(df.A.dg.assort, df.B.dg.assort, df.C.dg.assort, df.D.dg.assort, df.rest.dg.assort)
  
  return(df.fi)
}
