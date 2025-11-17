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
    df.merge.tmp <- as.data.frame(df.merge[i,1]) %>% dplyr::mutate(id.x = df.merge[i,4], isolate.y = df.merge[i,6], id.y = df.merge[i,9], sum = df.merge[i,2] + df.merge[i,7])
    names(df.merge.tmp) <- c("isolate.x", "id.x", "isolate.y", "id.y","sum")
    df.fi <- rbind(df.fi, df.merge.tmp)
  }
  return(df.fi)
}

Execution.kmer.count.matrix <- function(df) {
  df.output.1 <- kmer.count.matrix(df) %>% dplyr::mutate(subset = "1")
  df.output.2 <- kmer.count.matrix(df) %>% dplyr::mutate(subset = "2")
  df.output.3 <- kmer.count.matrix(df) %>% dplyr::mutate(subset = "3")
  df.output.4 <- kmer.count.matrix(df) %>% dplyr::mutate(subset = "4")
  df.output.5 <- kmer.count.matrix(df) %>% dplyr::mutate(subset = "5")
  df.output.6 <- kmer.count.matrix(df) %>% dplyr::mutate(subset = "6")
  df.output.7 <- kmer.count.matrix(df) %>% dplyr::mutate(subset = "7")
  df.output.8 <- kmer.count.matrix(df) %>% dplyr::mutate(subset = "8")
  df.output.9 <- kmer.count.matrix(df) %>% dplyr::mutate(subset = "9")
  df.output.10 <- kmer.count.matrix(df) %>% dplyr::mutate(subset = "10")
  
  df.output <- dplyr::bind_rows(df.output.1, df.output.2, df.output.3, df.output.4, df.output.5, df.output.6, df.output.7, df.output.8, df.output.9, df.output.10)
  
  return(df.output)
}

Execution.kmer.count.matrix.into.isolate <- function(df) {
  df.A <- df %>% dplyr::filter(subtype == "A")
  df.B <- df %>% dplyr::filter(subtype == "B")
  df.C <- df %>% dplyr::filter(subtype == "C")
  df.D <- df %>% dplyr::filter(subtype == "D")
  df.rest <- df %>% dplyr::filter(subtype == "rest")
  
  df.A.output <- Execution.kmer.count.matrix(df.A) %>% dplyr::mutate(subtype = "A")
  df.A.output$weight <- c(1-(df.A.output$sum)/sum(df.A.output$sum))
  
  df.B.output <- Execution.kmer.count.matrix(df.B) %>% dplyr::mutate(subtype = "B")
  df.B.output$weight <- c(1-(df.B.output$sum)/sum(df.B.output$sum))
  
  df.C.output <- Execution.kmer.count.matrix(df.C) %>% dplyr::mutate(subtype = "C")
  df.C.output$weight <- c(1-(df.C.output$sum)/sum(df.C.output$sum))
  
  df.D.output <- Execution.kmer.count.matrix(df.D) %>% dplyr::mutate(subtype = "D")
  df.D.output$weight <- c(1-(df.D.output$sum)/sum(df.D.output$sum))
  
  df.rest.output <- Execution.kmer.count.matrix(df.rest) %>% dplyr::mutate(subtype = "rest")
  df.rest.output$weight <- c(1-(df.rest.output$sum)/sum(df.rest.output$sum))
  
  df.subtyps.output <- dplyr::bind_rows(df.A.output, df.B.output, df.C.output, df.D.output, df.rest.output)
  
  
  
  return(df.subtyps.output)
}

plot.network.isolate <- function(df) {
  df.input <- df %>% dplyr::filter(subset == "1") %>% dplyr::select(isolate.x, isolate.y, weight)
  
  df.isolate.x <- df %>% dplyr::filter(subset == "1") %>% dplyr::select(isolate.x, subtype) %>% dplyr::rename(isolate = isolate.x)
  df.isolate.y <- df %>% dplyr::filter(subset == "1") %>% dplyr::select(isolate.y, subtype) %>% dplyr::rename(isolate = isolate.y)
  df.isolate.subtype.pair <- dplyr::bind_rows(df.isolate.x, df.isolate.y) %>% unique() 
  node.isolate <- df.isolate.subtype.pair %>% dplyr::select(isolate) %>% dplyr::mutate(count = 1)
  node.isolate.agg <- aggregate(count ~ isolate, node.isolate, sum) %>% arrange(desc(count))
  
  node.isolate.count.1 <- node.isolate.agg %>% dplyr::filter(count == 1)
  node.isolate.count.1.sub <- merge(node.isolate.count.1, df.isolate.subtype.pair, by = "isolate")
  
  node.isolate.count.23 <- node.isolate.agg %>% dplyr::filter(count != 1) %>% dplyr::mutate(subtype = case_when(count == 2 ~ "double",
                                                                                                                count == 3 ~ "triple"))
  
  node.fi <- bind_rows(node.isolate.count.1.sub, node.isolate.count.23)
  
  Richtung.subtype <- c("A", "B", "C", "D", "rest", "double", "triple")
  node.fi$subtype <- factor(node.fi$subtype, levels = Richtung.subtype)
  
  network.igraph <- graph_from_data_frame(d = df.input, vertices = node.fi, directed = F)
  
  color <- c("green4", "orange", "red3", "purple3", "navy", "grey", "grey4")
  color.subtype <- color[as.numeric(as.factor(V(network.igraph)$subtype))]
  
  plot(network.igraph, vertex.size = 10, vertex.label = NA, vertex.color = color.subtype)
}
