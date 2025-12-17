{r FUN::isolate graph matrix::mislabel}
# function
isolate.graph.matrix.mislabel <- function(df) {
  set.seed(42)
  
  total.count <- as.data.frame(rnorm(nrow(df), 0, 3))
  names(total.count) <- c("sum")
  total.count <- total.count %>% dplyr::mutate(series.id = c(1:nrow(df)))
  
  df$sum <- NULL
  df <- df %>% dplyr::mutate(series.id = c(1:nrow(df)))
  df.random <- dplyr::bind_cols(df, total.count, by = "series.id") %>% dplyr::select(isolate.x, id.x, isolate.y, id.y, sum, subset, subtype) 
  
  return(df.random)
}

Execution.kmer.count.matrix.into.isolate.mislabel <- function(df) {
  df.A <- df %>% dplyr::filter(subtype == "A")
  df.B <- df %>% dplyr::filter(subtype == "B")
  df.C <- df %>% dplyr::filter(subtype == "C")
  df.D <- df %>% dplyr::filter(subtype == "D")
  df.rest <- df %>% dplyr::filter(subtype == "rest")
  
  df.A.output <- isolate.graph.matrix.mislabel(df.A) %>% dplyr::mutate(subtype = "A")
  df.B.output <- isolate.graph.matrix.mislabel(df.B) %>% dplyr::mutate(subtype = "B")
  df.C.output <- isolate.graph.matrix.mislabel(df.C) %>% dplyr::mutate(subtype = "C")
  df.D.output <- isolate.graph.matrix.mislabel(df.D) %>% dplyr::mutate(subtype = "D")
  df.rest.output <- isolate.graph.matrix.mislabel(df.rest) %>% dplyr::mutate(subtype = "rest")
  
  df.subtypes.output <- dplyr::bind_rows(df.A.output, df.B.output, df.C.output, df.D.output, df.rest.output)
  
  return(df.subtypes.output)
}

# execution
subsets <- c("1","2", "3", "4", "5", "6", "7", "8", "9", "10")

isolate.kmer.count.mx.DNA.13.mislabel.list <- lapply(subsets, function(i) {
  Execution.kmer.count.matrix.into.isolate.mislabel(isolate.kmer.count.mx.DNA.13) %>% dplyr::mutate(subset = i)
})

isolate.kmer.count.mx.DNA.15.mislabel.list <- lapply(subsets, function(i) {
  Execution.kmer.count.matrix.into.isolate.mislabel(isolate.kmer.count.mx.DNA.15) %>% dplyr::mutate(subset = i)
})

isolate.kmer.count.mx.RNA.13.mislabel.list <- lapply(subsets, function(i) {
  Execution.kmer.count.matrix.into.isolate.mislabel(isolate.kmer.count.mx.RNA.13) %>% dplyr::mutate(subset = i)
})

isolate.kmer.count.mx.RNA.15.mislabel.list <- lapply(subsets, function(i) {
  Execution.kmer.count.matrix.into.isolate.mislabel(isolate.kmer.count.mx.RNA.15) %>% dplyr::mutate(subset = i)
})

{r reshape graph matrix::mislabel}
# function
reshape.graph.matrix.isolate <- function(df) {
  df_ <- df %>% dplyr::select(isolate.x, isolate.y, sum) %>% filter(sum != 0)
  df.complete <- df_ %>% complete(isolate.x, isolate.y, fill = list(sum = 0))
  df.wide <- dcast(df.complete, isolate.x ~ isolate.y, value.var = "sum", fun.aggregate = sum)
  rownames(df.wide) <- df.wide$isolate.x
  df.wide$isolate.x <- NULL
  output.mx <- as.matrix(df.wide)
  
  return(output.mx)
}

# execution
graph.matrix.DNA.13.mislabel.isolate.list <- lapply(1:10, function(i) {
  reshape.graph.matrix.isolate(isolate.kmer.count.mx.DNA.13.mislabel.list[[i]])
})

graph.matrix.DNA.15.mislabel.isolate.list <- lapply(1:10, function(i) {
  reshape.graph.matrix.isolate(isolate.kmer.count.mx.DNA.15.mislabel.list[[i]])
})

graph.matrix.RNA.13.mislabel.isolate.list <- lapply(1:10, function(i) {
  reshape.graph.matrix.isolate(isolate.kmer.count.mx.RNA.13.mislabel.list[[i]])
})

graph.matrix.RNA.15.mislabel.isolate.list <- lapply(1:10, function(i) {
  reshape.graph.matrix.isolate(isolate.kmer.count.mx.RNA.15.mislabel.list[[i]])
})

{r prepare matrix coordinate for the row::mislabel}
# function
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

# execution
mx.coordinate.DNA.13.mislabel.isolate.row.list <- lapply(1:10, function(i) {
  matrix.coordinate.isolate(isolate.kmer.count.mx.DNA.13.list[[i]], graph.matrix.DNA.13.mislabel.isolate.list[[i]])
})

mx.coordinate.DNA.15.mislabel.isolate.row.list <- lapply(1:10, function(i) {
  matrix.coordinate.isolate(isolate.kmer.count.mx.DNA.15.list[[i]], graph.matrix.DNA.15.mislabel.isolate.list[[i]])
})

mx.coordinate.RNA.13.mislabel.isolate.row.list <- lapply(1:10, function(i) {
  matrix.coordinate.isolate(isolate.kmer.count.mx.RNA.13.list[[i]], graph.matrix.RNA.13.mislabel.isolate.list[[i]])
})

mx.coordinate.RNA.15.mislabel.isolate.row.list <- lapply(1:10, function(i) {
  matrix.coordinate.isolate(isolate.kmer.count.mx.RNA.15.list[[i]], graph.matrix.RNA.15.mislabel.isolate.list[[i]])
})

{r Random walk}
# function
simulate_random_walk.subtype <- function(graph_matrix, start_node, num_steps) {
  
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
