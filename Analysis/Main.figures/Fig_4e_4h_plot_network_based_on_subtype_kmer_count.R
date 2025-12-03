plot.network <- function(df) {
  df.input <- df %>% dplyr::filter(subset == "1") %>% dplyr::select(kmer.x, kmer.y, distance)
  
  df.kmer.x <- df %>% dplyr::filter(subset == "1") %>% dplyr::select(kmer.x, subtype) %>% dplyr::rename(kmer = kmer.x)
  df.kmer.y <- df %>% dplyr::filter(subset == "1") %>% dplyr::select(kmer.y, subtype) %>% dplyr::rename(kmer = kmer.y)
  df.kmer.subtype.pair <- dplyr::bind_rows(df.kmer.x, df.kmer.y) %>% unique() 
  node.kmer <- df.kmer.subtype.pair %>% dplyr::select(kmer) %>% dplyr::mutate(count = 1)
  node.kmer.agg <- aggregate(count ~ kmer, node.kmer, sum) %>% arrange(desc(count))
  
  node.kmer.count.1 <- node.kmer.agg %>% dplyr::filter(count == 1)
  node.kmer.count.1.sub <- merge(node.kmer.count.1, df.kmer.subtype.pair, by = "kmer")
  
  node.kmer.count.23 <- node.kmer.agg %>% dplyr::filter(count != 1) %>% dplyr::mutate(subtype = case_when(count == 2 ~ "double",
                                                                                                          count == 3 ~ "triple"))
  
  node.fi <- bind_rows(node.kmer.count.1.sub, node.kmer.count.23)
  
  Richtung.subtype <- c("A", "B", "C", "D", "rest", "double", "triple")
  node.fi$subtype <- factor(node.fi$subtype, levels = Richtung.subtype)
  
  network.igraph <- graph_from_data_frame(d = df.input, vertices = node.fi, directed = F)
  
  color <- c("green4", "orange", "red3", "purple3", "navy", "lightgrey", "grey3")
  color.subtype <- color[as.numeric(as.factor(V(network.igraph)$subtype))]
  
  plot(network.igraph, vertex.size = 5, vertex.label = NA, vertex.color = color.subtype)
}
