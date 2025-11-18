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

Execution.kmer.count.matrix.into.isolate.weight <- function(df) {
  df.A <- df %>% dplyr::filter(subtype == "A")
  df.B <- df %>% dplyr::filter(subtype == "B")
  df.C <- df %>% dplyr::filter(subtype == "C")
  df.D <- df %>% dplyr::filter(subtype == "D")
  df.rest <- df %>% dplyr::filter(subtype == "rest")
  
  df.A.output <- Execution.kmer.count.matrix(df.A) %>% dplyr::mutate(subtype = "A")
  df.A.output$distance <- c((df.A.output$sum)/sum(df.A.output$sum))
  
  df.B.output <- Execution.kmer.count.matrix(df.B) %>% dplyr::mutate(subtype = "B")
  df.B.output$distance <- c((df.B.output$sum)/sum(df.B.output$sum))
  
  df.C.output <- Execution.kmer.count.matrix(df.C) %>% dplyr::mutate(subtype = "C")
  df.C.output$distance <- c((df.C.output$sum)/sum(df.C.output$sum))
  
  df.D.output <- Execution.kmer.count.matrix(df.D) %>% dplyr::mutate(subtype = "D")
  df.D.output$distance <- c((df.D.output$sum)/sum(df.D.output$sum))
  
  df.rest.output <- Execution.kmer.count.matrix(df.rest) %>% dplyr::mutate(subtype = "rest")
  df.rest.output$distance <- c((df.rest.output$sum)/sum(df.rest.output$sum))
  
  df.subtyps.output <- dplyr::bind_rows(df.A.output, df.B.output, df.C.output, df.D.output, df.rest.output)
  
  
  
  return(df.subtyps.output)
}

Euclidean.heatmap.plot <- function(df, df.kmer) {
  df <- df %>% dplyr::filter(subset == "1") %>% dplyr::select(isolate.x, distance) %>% dplyr::rename(isolate = isolate.x)
  
  df.kmer.merge <- merge(df, df.kmer, by = "isolate") %>% dplyr::select(isolate, group, distance, kmer.count)
  
  df.kmer.merge.agg <- aggregate(. ~ isolate + group, df.kmer.merge, mean) %>% arrange(desc(isolate))
  
  df.kmer.merge.agg.eucl <- dist(x = df.kmer.merge.agg[,3:4], method = "euclidean")
  
  df.kmer.merge.agg.eucl.mx <- data.matrix(df.kmer.merge.agg.eucl)
  
  col_fun <- colorRamp2(c(0, 2.5, 5), c("yellow2", "black", "midnightblue"))
  
  df.anno <- df.kmer.merge.agg %>% dplyr::select(isolate, group)
  rownames(df.anno) <- df.anno$isolate
  df.anno$isolate <- NULL
  
  ha.anno.row <- rowAnnotation(df = df.anno, col = list("group" = c("A" = "green4", "B" = "orange", "C" = "red3", "D" = "purple3", "rest" = "navy")))
  
  ha.anno.col <- HeatmapAnnotation(df = df.anno, col = list("group" = c("A" = "green4", "B" = "orange", "C" = "red3", "D" = "purple3", "rest" = "navy")))
  
  Heatmap(df.kmer.merge.agg.eucl.mx, col = col_fun, clustering_distance_columns = "euclidean", clustering_distance_rows = "euclidean", show_row_names = F, show_column_names = F, right_annotation = ha.anno.row, bottom_annotation = ha.anno.col)
}
