{r FUN::subtype graph matrix::mislabel}
# function
subtype.graph.matrix.mislabel <- function(df) {
  set.seed(42)
  
  total.count <- as.data.frame(rnorm(nrow(df), 0, 3))
  names(total.count) <- c("sum")
  total.count <- total.count %>% dplyr::mutate(series.id = c(1:nrow(df)))
  
  df$sum <- NULL
  df <- df %>% dplyr::mutate(series.id = c(1:nrow(df)))
  df.random <- dplyr::bind_cols(df, total.count, by = "series.id") %>% dplyr::select(kmer.x, id.x, kmer.y, id.y, sum, subset, subtype) 
  
  return(df.random)
}

Execution.kmer.count.matrix.into.subtype.mislabel <- function(df) {
  df.A <- df %>% dplyr::filter(subtype == "A")
  df.B <- df %>% dplyr::filter(subtype == "B")
  df.C <- df %>% dplyr::filter(subtype == "C")
  df.D <- df %>% dplyr::filter(subtype == "D")
  df.rest <- df %>% dplyr::filter(subtype == "rest")
  
  df.A.output <- subtype.graph.matrix.mislabel(df.A) %>% dplyr::mutate(subtype = "A")
  df.B.output <- subtype.graph.matrix.mislabel(df.B) %>% dplyr::mutate(subtype = "B")
  df.C.output <- subtype.graph.matrix.mislabel(df.C) %>% dplyr::mutate(subtype = "C")
  df.D.output <- subtype.graph.matrix.mislabel(df.D) %>% dplyr::mutate(subtype = "D")
  df.rest.output <- subtype.graph.matrix.mislabel(df.rest) %>% dplyr::mutate(subtype = "rest")
  
  df.subtypes.output <- dplyr::bind_rows(df.A.output, df.B.output, df.C.output, df.D.output, df.rest.output)
  
  return(df.subtypes.output)
}

# execution
subsets <- c("1","2", "3", "4", "5", "6", "7", "8", "9", "10")

subtype.kmer.count.mx.DNA.13.mislabel.list <- lapply(subsets, function(i) {
  Execution.kmer.count.matrix.into.subtype.mislabel(subtype.kmer.count.mx.DNA.13) %>% dplyr::mutate(subset = i)
})

subtype.kmer.count.mx.DNA.15.mislabel.list <- lapply(subsets, function(i) {
  Execution.kmer.count.matrix.into.subtype.mislabel(subtype.kmer.count.mx.DNA.15) %>% dplyr::mutate(subset = i)
})

subtype.kmer.count.mx.RNA.13.mislabel.list <- lapply(subsets, function(i) {
  Execution.kmer.count.matrix.into.subtype.mislabel(subtype.kmer.count.mx.RNA.13) %>% dplyr::mutate(subset = i)
})

subtype.kmer.count.mx.RNA.15.mislabel.list <- lapply(subsets, function(i) {
  Execution.kmer.count.matrix.into.subtype.mislabel(subtype.kmer.count.mx.RNA.15) %>% dplyr::mutate(subset = i)
})

{r reshape graph matrix::mislabel}
#function
reshape.graph.matrix <- function(df) {
  df_ <- df %>% dplyr::select(id.x, id.y, sum) %>% filter(sum != 0)
  df.complete <- df_ %>% complete(id.x, id.y, fill = list(sum = 0))
  df.wide <- dcast(df.complete, id.x ~ id.y, value.var = "sum", fun.aggregate = sum)
  rownames(df.wide) <- df.wide$id.x
  df.wide$id.x <- NULL
  output.mx <- as.matrix(df.wide)
  
  return(output.mx)
}

# execution
graph.matrix.DNA.13.mislabel.list <- lapply(1:10, function(i) {
  reshape.graph.matrix(subtype.kmer.count.mx.DNA.13.mislabel.list[[i]])
})

graph.matrix.DNA.15.mislabel.list <- lapply(1:10, function(i) {
  reshape.graph.matrix(subtype.kmer.count.mx.DNA.15.mislabel.list[[i]])
})

graph.matrix.RNA.13.mislabel.list <- lapply(1:10, function(i) {
  reshape.graph.matrix(subtype.kmer.count.mx.RNA.13.mislabel.list[[i]])
})

graph.matrix.RNA.15.mislabel.list <- lapply(1:10, function(i) {
  reshape.graph.matrix(subtype.kmer.count.mx.RNA.15.mislabel.list[[i]])
})

graph.matrix.DNA.13.mislabel.list <- lapply(1:10, function(i) {
  reshape.graph.matrix(subtype.kmer.count.mx.DNA.13.mislabel.list[[i]])
})

graph.matrix.DNA.15.mislabel.list <- lapply(1:10, function(i) {
  reshape.graph.matrix(subtype.kmer.count.mx.DNA.15.mislabel.list[[i]])
})

graph.matrix.RNA.13.mislabel.list <- lapply(1:10, function(i) {
  reshape.graph.matrix(subtype.kmer.count.mx.RNA.13.mislabel.list[[i]])
})

graph.matrix.RNA.15.mislabel.list <- lapply(1:10, function(i) {
  reshape.graph.matrix(subtype.kmer.count.mx.RNA.15.mislabel.list[[i]])
})

{r prepare matrix coordinate for the row::mislabel}
#function
node.list <- function(mx) {
  mx.x <- mx %>% dplyr::select(id.x) %>% dplyr::rename(node = id.x)
  mx.y <- mx %>% dplyr::select(id.y) %>% dplyr::rename(node = id.y)
  
  node.list <- dplyr::bind_rows(mx.x, mx.y) %>% unique()
  node.list$node <- as.character(node.list$node)
  return(node.list)
}

matrix.coordinate <- function(mx, graph.mx) {
  df.row <- data.frame(dimnames(graph.mx)[1])
  colnames(df.row) <- c("node")
  df.row <- dplyr::mutate(df.row, order = "row", id = c(1:nrow(df.row)))
  
  mx.A <- mx %>% dplyr::filter(subtype == "A")
  node.A <- node.list(mx.A)
  row.A <- dplyr::inner_join(df.row, node.A, by = "node") %>% mutate(subtype = "A")
  
  mx.B <- mx %>% dplyr::filter(subtype == "B") 
  node.B <- node.list(mx.B)
  row.B <- dplyr::inner_join(df.row, node.B, by = "node") %>% mutate(subtype = "B")
  
  mx.C <- mx %>% dplyr::filter(subtype == "C")
  node.C <- node.list(mx.C)
  row.C <- dplyr::inner_join(df.row, node.C, by = "node") %>% mutate(subtype = "C")
  
  mx.D <- mx %>% dplyr::filter(subtype == "D")
  node.D <- node.list(mx.D)
  row.D <- dplyr::inner_join(df.row, node.D, by = "node") %>% mutate(subtype = "D") 
  
  mx.R <- mx %>% dplyr::filter(subtype == "rest")
  node.R <- node.list(mx.R)
  row.R <- dplyr::inner_join(df.row, node.R, by = "node") %>% mutate(subtype = "R")
  
  df.complete <- dplyr::bind_rows(row.A, row.B, row.C, row.D, row.R)
  
  return(df.complete)
}

# execution
mx.coordinate.DNA.13.mislabel.row.list <- lapply(1:10, function(i) {
  matrix.coordinate(subtype.kmer.count.mx.DNA.13.list[[i]], graph.matrix.DNA.13.mislabel.list[[i]])
})

mx.coordinate.DNA.15.mislabel.row.list <- lapply(1:10, function(i) {
  matrix.coordinate(subtype.kmer.count.mx.DNA.15.list[[i]], graph.matrix.DNA.15.mislabel.list[[i]])
})

mx.coordinate.RNA.13.mislabel.row.list <- lapply(1:10, function(i) {
  matrix.coordinate(subtype.kmer.count.mx.RNA.13.list[[i]], graph.matrix.RNA.13.mislabel.list[[i]])
})

mx.coordinate.RNA.15.mislabel.row.list <- lapply(1:10, function(i) {
  matrix.coordinate(subtype.kmer.count.mx.RNA.15.list[[i]], graph.matrix.RNA.15.mislabel.list[[i]])
})

{r Random walk}
#function
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
