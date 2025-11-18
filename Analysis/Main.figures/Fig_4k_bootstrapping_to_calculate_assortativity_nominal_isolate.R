calculation.mismatch.assortativity.nominal <- function(df) {
  set.seed(42)
  
  total.count <- as.data.frame(sample(1:20000, nrow(df), replace =T))
  names(total.count) <- c("total.count")
  
  df.random <- dplyr::bind_cols(df, total.count) %>% dplyr::select(isolate.x, id.x, isolate.y, id.y, total.count, subset, subtype) 
  df.random <- df.random %>% dplyr::mutate(distance = 1-(df.random$total.count)/sum(df.random$total.count))
  
  df.input <- df.random %>% dplyr::select(isolate.x, isolate.y, distance)
  
  df.isolate.x <- df.random %>% dplyr::select(isolate.x, subtype) %>% dplyr::rename(isolate = isolate.x) %>% unique()
  df.isolate.y <- df.random %>% dplyr::select(isolate.y, subtype) %>% dplyr::rename(isolate = isolate.y) %>% unique()
  df.isolate.subtype.pair <- dplyr::bind_rows(df.isolate.x, df.isolate.y) 
  node.isolate <- df.isolate.subtype.pair %>% dplyr::select(isolate) %>% dplyr::mutate(count = 1)
  node.isolate.agg <- aggregate(count ~ isolate, node.isolate, sum) %>% arrange(desc(count))
  
  node.isolate.count.1 <- node.isolate.agg %>% dplyr::filter(count == 1)
  node.isolate.count.1.sub <- merge(node.isolate.count.1, df.isolate.subtype.pair, by = "isolate")
  
  node.isolate.count.23 <- node.isolate.agg %>% dplyr::filter(count != 1) %>% dplyr::mutate(subtype = case_when(count == 2 ~ "double",
                                                                                                                count == 3 | count > 3 ~ "triple"))
  
  node.fi <- bind_rows(node.isolate.count.1.sub, node.isolate.count.23)
  
  Richtung.subtype <- c("A", "B", "C", "D", "rest", "double", "triple")
  node.fi$subtype <- factor(node.fi$subtype, levels = Richtung.subtype)
  
  network.igraph <- graph_from_data_frame(d = df.input, vertices = node.fi, directed = F)
  assort.nominal <- assortativity_nominal(network.igraph, as.integer(as.factor(V(network.igraph)$subtype)))
  
  return(assort.nominal)
}

calculation.mismatch.assortativity.nominal.subset <- function(df) {
  
  df$distance <- NULL
  df$sum <- NULL
  
  df.subset.1 <- df %>% dplyr::filter(subset == "1")
  df.subset.1.assort.nominal <- calculation.mismatch.assortativity.nominal(df.subset.1)
  
  df.subset.2 <- df %>% dplyr::filter(subset == "2")
  df.subset.2.assort.nominal <- calculation.mismatch.assortativity.nominal(df.subset.2)
  
  df.subset.3 <- df %>% dplyr::filter(subset == "3")
  df.subset.3.assort.nominal <- calculation.mismatch.assortativity.nominal(df.subset.3)
  
  df.subset.4 <- df %>% dplyr::filter(subset == "4")
  df.subset.4.assort.nominal <- calculation.mismatch.assortativity.nominal(df.subset.4)
  
  df.subset.5 <- df %>% dplyr::filter(subset == "5")
  df.subset.5.assort.nominal <- calculation.mismatch.assortativity.nominal(df.subset.5)
  
  df.subset.6 <- df %>% dplyr::filter(subset == "6")
  df.subset.6.assort.nominal <- calculation.mismatch.assortativity.nominal(df.subset.6)
  
  df.subset.7 <- df %>% dplyr::filter(subset == "7")
  df.subset.7.assort.nominal <- calculation.mismatch.assortativity.nominal(df.subset.7)
  
  df.subset.8 <- df %>% dplyr::filter(subset == "8")
  df.subset.8.assort.nominal <- calculation.mismatch.assortativity.nominal(df.subset.8)
  
  df.subset.9 <- df %>% dplyr::filter(subset == "9")
  df.subset.9.assort.nominal <- calculation.mismatch.assortativity.nominal(df.subset.9)
  
  df.subset.10 <- df %>% dplyr::filter(subset == "10")
  df.subset.10.assort.nominal <- calculation.mismatch.assortativity.nominal(df.subset.10)
  
  value <- c(df.subset.1.assort.nominal, df.subset.2.assort.nominal, df.subset.3.assort.nominal, df.subset.4.assort.nominal, df.subset.5.assort.nominal, df.subset.6.assort.nominal, df.subset.7.assort.nominal, df.subset.8.assort.nominal, df.subset.9.assort.nominal, df.subset.10.assort.nominal)
  subset <- c(1:10)
  sample <- rep(c("Mismatch"), 10)
  
  df.subset.output <- data.frame(value, subset, sample)
  df.subset.output$subset <- as.character(df.subset.output$subset)
  
  return(df.subset.output)
}
