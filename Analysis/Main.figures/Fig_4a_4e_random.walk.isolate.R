```{r FUN::isolate graph matrix}
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

reshape.graph.matrix.isolate <- function(df) {
  df_ <- df %>% dplyr::select(isolate.x, isolate.y, sum) %>% filter(sum != 0)
  df.complete <- df_ %>% complete(isolate.x, isolate.y, fill = list(sum = 0))
  df.wide <- dcast(df.complete, isolate.x ~ isolate.y, value.var = "sum", fun.aggregate = sum)
  rownames(df.wide) <- df.wide$isolate.x
  df.wide$isolate.x <- NULL
  output.mx <- as.matrix(df.wide)
  
  return(output.mx)
}

isolate.list <- function(mx) {
  mx.x <- mx %>% dplyr::select(isolate.x) %>% dplyr::rename(isolate = isolate.x)
  mx.y <- mx %>% dplyr::select(isolate.y) %>% dplyr::rename(isolate = isolate.y)
  
  isolate.list <- dplyr::bind_rows(mx.x, mx.y) %>% unique()
  isolate.list$isolate <- as.character(isolate.list$isolate)
  return(isolate.list)
}

matrix.coordinate.isolate <- function(mx, graph.mx) {
  df.row <- data.frame(dimnames(graph.mx)[1])
  colnames(df.row) <- c("isolate")
  df.row <- dplyr::mutate(df.row, order = "row", id = c(1:nrow(df.row)))
  
  mx.A <- mx %>% dplyr::filter(subtype == "A")
  isolate.A <- isolate.list(mx.A)
  row.A <- dplyr::inner_join(df.row, isolate.A, by = "isolate") %>% mutate(subtype = "A")
  
  mx.B <- mx %>% dplyr::filter(subtype == "B") 
  isolate.B <- isolate.list(mx.B)
  row.B <- dplyr::inner_join(df.row, isolate.B, by = "isolate") %>% mutate(subtype = "B")
  
  mx.C <- mx %>% dplyr::filter(subtype == "C")
  isolate.C <- isolate.list(mx.C)
  row.C <- dplyr::inner_join(df.row, isolate.C, by = "isolate") %>% mutate(subtype = "C")
  
  mx.D <- mx %>% dplyr::filter(subtype == "D")
  isolate.D <- isolate.list(mx.D)
  row.D <- dplyr::inner_join(df.row, isolate.D, by = "isolate") %>% mutate(subtype = "D") 
  
  mx.R <- mx %>% dplyr::filter(subtype == "rest")
  isolate.R <- isolate.list(mx.R)
  row.R <- dplyr::inner_join(df.row, isolate.R, by = "isolate") %>% mutate(subtype = "R")
  
  df.complete <- dplyr::bind_rows(row.A, row.B, row.C, row.D, row.R)
  
  return(df.complete)
}

mx.coordinate.DNA.13.isolate.row.list <- lapply(1:10, function(i) {
  matrix.coordinate.isolate(isolate.kmer.count.mx.DNA.13.list[[i]], graph.matrix.DNA.13.isolate.list[[i]])
})

mx.coordinate.DNA.15.isolate.row.list <- lapply(1:10, function(i) {
  matrix.coordinate.isolate(isolate.kmer.count.mx.DNA.15.list[[i]], graph.matrix.DNA.15.isolate.list[[i]])
})

mx.coordinate.RNA.13.isolate.row.list <- lapply(1:10, function(i) {
  matrix.coordinate.isolate(isolate.kmer.count.mx.RNA.13.list[[i]], graph.matrix.RNA.13.isolate.list[[i]])
})

mx.coordinate.RNA.15.isolate.row.list <- lapply(1:10, function(i) {
  matrix.coordinate.isolate(isolate.kmer.count.mx.RNA.15.list[[i]], graph.matrix.RNA.15.isolate.list[[i]])
})

simulate_random_walk.subtype <- function(graph_matrix, start_node, num_steps) {
  
  #graph_matrix[graph_matrix != 0] <- 1
  #graph_matrix.t <- t(graph_matrix)
  
  #num_nodes <- nrow(graph_matrix)
  path <- numeric(num_steps)
  path[1] <- start_node
  current_node <- start_node
  
  for (i in 1:num_steps) {
    # Find neighbors of the current node (row)
    neighbors <- which(graph_matrix[current_node, ] != 0)
    
    # If no neighbors, stop the walk (or define another rule)
    if (length(neighbors) == 0) {
      break
    }
    
    # Randomly select the next node from neighbors
    current_node <- sample(neighbors, 1)
    path[i + 1] <- current_node
  }
  return(path)
}
