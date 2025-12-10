{r subtype graph matrix}
# function
subtype.graph.matrix <- function(df) {
  df$id <- c(1:nrow(df)) # indexing K-mers
  df.samp <- sample_n(df, 100, replace = T)
  df1 <- df.samp %>% dplyr::select(kmer, kmer.count, subtype, id) %>% dplyr::mutate(temp = "1")
  df2 <- df.samp %>% dplyr::select(kmer, kmer.count, subtype, id) %>% dplyr::mutate(temp = "1")
  
  df.merge <- dplyr::full_join(df1, df2, by = "temp", relationship = "many-to-many")
  
  df.fi <- data.frame()
  
  for(i in 1:nrow(df.merge)) {
    if (isTRUE(df.merge[i,1] == df.merge[i,6])) {
      next
    }
    df.merge.tmp <- as.data.frame(df.merge[i,1]) %>% dplyr::mutate(id.x = df.merge[i,4], kmer.y = df.merge[i,6], id.y = df.merge[i,9], sum = df.merge[i,2] + df.merge[i,7])
    names(df.merge.tmp) <- c("kmer.x", "id.x", "kmer.y", "id.y","sum")
    df.fi <- rbind(df.fi, df.merge.tmp)
  }
  return(df.fi)
}

Execution.kmer.count.matrix.into.subtype <- function(df) {
  df.A <- df %>% dplyr::filter(subtype == "A")
  df.B <- df %>% dplyr::filter(subtype == "B")
  df.C <- df %>% dplyr::filter(subtype == "C")
  df.D <- df %>% dplyr::filter(subtype == "D")
  df.rest <- df %>% dplyr::filter(subtype == "rest")
  
  df.A.output <- subtype.graph.matrix(df.A) %>% dplyr::mutate(subtype = "A")
  df.B.output <- subtype.graph.matrix(df.B) %>% dplyr::mutate(subtype = "B")
  df.C.output <- subtype.graph.matrix(df.C) %>% dplyr::mutate(subtype = "C")
  df.D.output <- subtype.graph.matrix(df.D) %>% dplyr::mutate(subtype = "D")
  df.rest.output <- subtype.graph.matrix(df.rest) %>% dplyr::mutate(subtype = "rest")
  
  df.subtypes.output <- dplyr::bind_rows(df.A.output, df.B.output, df.C.output, df.D.output, df.rest.output)
  
  return(df.subtypes.output)
}

# execution
subsets <- c("1","2", "3", "4", "5", "6", "7", "8", "9", "10")

subtype.kmer.count.mx.DNA.13.list <- lapply(subsets, function(i) {
  Execution.kmer.count.matrix.into.subtype(subtype.kmer.count.DNA.13) %>% dplyr::mutate(subset = i)
})

subtype.kmer.count.mx.DNA.15.list <- lapply(subsets, function(i) {
  Execution.kmer.count.matrix.into.subtype(subtype.kmer.count.DNA.15) %>% dplyr::mutate(subset = i)
})

subtype.kmer.count.mx.RNA.13.list <- lapply(subsets, function(i) {
  Execution.kmer.count.matrix.into.subtype(subtype.kmer.count.RNA.13) %>% dplyr::mutate(subset = i)
})

subtype.kmer.count.mx.RNA.15.list <- lapply(subsets, function(i) {
  Execution.kmer.count.matrix.into.subtype(subtype.kmer.count.RNA.15) %>% dplyr::mutate(subset = i)
})

{r reshape graph matrix}
# function
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
graph.matrix.DNA.13.list <- lapply(1:10, function(i) {
  reshape.graph.matrix(subtype.kmer.count.mx.DNA.13.list[[i]])
})

graph.matrix.DNA.15.list <- lapply(1:10, function(i) {
  reshape.graph.matrix(subtype.kmer.count.mx.DNA.15.list[[i]])
})

graph.matrix.RNA.13.list <- lapply(1:10, function(i) {
  reshape.graph.matrix(subtype.kmer.count.mx.RNA.13.list[[i]])
})

graph.matrix.RNA.15.list <- lapply(1:10, function(i) {
  reshape.graph.matrix(subtype.kmer.count.mx.RNA.15.list[[i]])
})

{r prepare matrix coordinate, row direction}
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
mx.coordinate.DNA.13.row.list <- lapply(1:10, function(i) {
  matrix.coordinate(subtype.kmer.count.mx.DNA.13.list[[i]], graph.matrix.DNA.13.list[[i]])
})

mx.coordinate.DNA.15.row.list <- lapply(1:10, function(i) {
  matrix.coordinate(subtype.kmer.count.mx.DNA.15.list[[i]], graph.matrix.DNA.15.list[[i]])
})

mx.coordinate.RNA.13.row.list <- lapply(1:10, function(i) {
  matrix.coordinate(subtype.kmer.count.mx.RNA.13.list[[i]], graph.matrix.RNA.13.list[[i]])
})

mx.coordinate.RNA.15.row.list <- lapply(1:10, function(i) {
  matrix.coordinate(subtype.kmer.count.mx.RNA.15.list[[i]], graph.matrix.RNA.15.list[[i]])
})

{r Random walk}
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

walk.path.list <- function(graph.list) {
  start_node <- 2
  num_steps <- 10000
  walk.path.1 <- simulate_random_walk.subtype(graph.list[[1]], start_node, num_steps)
  
  start_node <- 12
  num_steps <- 10000
  walk.path.2 <- simulate_random_walk.subtype(graph.list[[2]], start_node, num_steps)
  
  start_node <- 12
  num_steps <- 10000
  walk.path.3 <- simulate_random_walk.subtype(graph.list[[3]], start_node, num_steps)
  
  start_node <- 9
  num_steps <- 10000
  walk.path.4 <- simulate_random_walk.subtype(graph.list[[4]], start_node, num_steps)
  
  start_node <- 6
  num_steps <- 10000
  walk.path.5 <- simulate_random_walk.subtype(graph.list[[5]], start_node, num_steps)
  
  start_node <- 1
  num_steps <- 10000
  walk.path.6 <- simulate_random_walk.subtype(graph.list[[6]], start_node, num_steps)
  
  start_node <- 3
  num_steps <- 10000
  walk.path.7 <- simulate_random_walk.subtype(graph.list[[7]], start_node, num_steps)
  
  start_node <- 13
  num_steps <- 10000
  walk.path.8 <- simulate_random_walk.subtype(graph.list[[8]], start_node, num_steps)
  
  start_node <- 4
  num_steps <- 10000
  walk.path.9 <- simulate_random_walk.subtype(graph.list[[9]], start_node, num_steps)
  
  start_node <- 7
  num_steps <- 10000
  walk.path.10 <- simulate_random_walk.subtype(graph.list[[10]], start_node, num_steps)
  
  walk.path.list <- list(walk.path.1, walk.path.2, walk.path.3, walk.path.4, walk.path.5, walk.path.6, walk.path.7, walk.path.8, walk.path.9, walk.path.10)
  
  return(walk.path.list)
}

# execution
walk_path.DNA.13.list <- walk.path.list(graph.matrix.DNA.13.list)
walk_path.DNA.15.list <- walk.path.list(graph.matrix.DNA.15.list)
walk_path.RNA.13.list <- walk.path.list(graph.matrix.RNA.13.list)
walk_path.RNA.15.list <- walk.path.list(graph.matrix.RNA.15.list)

{r Probability random walk}
# function
random.walk.prob <- function(path, coordinate) {
  id.A <- coordinate %>% dplyr::filter(subtype == "A")
  id.B <- coordinate %>% dplyr::filter(subtype == "B")
  id.C <- coordinate %>% dplyr::filter(subtype == "C")
  id.D <- coordinate %>% dplyr::filter(subtype == "D")
  id.R <- coordinate %>% dplyr::filter(subtype == "R")
  
  df.path <- as.data.frame(path)
  names(df.path) <- c("id")
  
  prob.A <- dim(dplyr::inner_join(df.path, id.A))[1]/dim(df.path)[1]
  prob.B <- dim(dplyr::inner_join(df.path, id.B))[1]/dim(df.path)[1]
  prob.C <- dim(dplyr::inner_join(df.path, id.C))[1]/dim(df.path)[1]
  prob.D <- dim(dplyr::inner_join(df.path, id.D))[1]/dim(df.path)[1]
  prob.R <- dim(dplyr::inner_join(df.path, id.R))[1]/dim(df.path)[1]
  
  df.prob <- data.frame(prob = c(prob.A, prob.B, prob.C, prob.D, prob.R), subtype = c("A", "B", "C", "D", "R"))
  
  return(df.prob)
}

combine.dataframe.list <- function(df) {
  df.split <- lapply(df, as.data.frame.list)
  df.merge <- do.call(what = "rbind", df.split)
  df.merge[is.na(df.merge)] <- 0
  
  return(df.merge)
}

# execution
## initiation A
df.prob.DNA.13.initiate.A <- combine.dataframe.list(lapply(1:10, function(i) {
  random.walk.prob(walk_path.DNA.13.list[[i]], mx.coordinate.DNA.13.row.list[[i]])
})) %>% dplyr::mutate(kmer.type = "DNA", kmer.length = 13)

df.prob.DNA.15.initiate.A <- combine.dataframe.list(lapply(1:10, function(i) {
  random.walk.prob(walk_path.DNA.15.list[[i]], mx.coordinate.DNA.15.row.list[[i]])
})) %>% dplyr::mutate(kmer.type = "DNA", kmer.length = 15)

df.prob.RNA.13.initiate.A <- combine.dataframe.list(lapply(1:10, function(i) {
  random.walk.prob(walk_path.RNA.13.list[[i]], mx.coordinate.RNA.13.row.list[[i]])
})) %>% dplyr::mutate(kmer.type = "RNA", kmer.length = 13)

df.prob.RNA.15.initiate.A <- combine.dataframe.list(lapply(1:10, function(i) {
  random.walk.prob(walk_path.RNA.15.list[[i]], mx.coordinate.RNA.15.row.list[[i]])
})) %>% dplyr::mutate(kmer.type = "RNA", kmer.length = 15)

df.prob.initiate.A <- bind_rows(df.prob.DNA.13.initiate.A, df.prob.DNA.15.initiate.A, df.prob.RNA.13.initiate.A, df.prob.RNA.15.initiate.A)

# initiation B
df.prob.DNA.13.initiate.B <- combine.dataframe.list(lapply(1:10, function(i) {
  random.walk.prob(walk_path.DNA.13.list[[i]], mx.coordinate.DNA.13.row.list[[i]])
})) %>% dplyr::mutate(kmer.type = "DNA", kmer.length = 13)

df.prob.DNA.15.initiate.B <- combine.dataframe.list(lapply(1:10, function(i) {
  random.walk.prob(walk_path.DNA.15.list[[i]], mx.coordinate.DNA.15.row.list[[i]])
})) %>% dplyr::mutate(kmer.type = "DNA", kmer.length = 15)

df.prob.RNA.13.initiate.B <- combine.dataframe.list(lapply(1:10, function(i) {
  random.walk.prob(walk_path.RNA.13.list[[i]], mx.coordinate.RNA.13.row.list[[i]])
})) %>% dplyr::mutate(kmer.type = "RNA", kmer.length = 13)

df.prob.RNA.15.initiate.B <- combine.dataframe.list(lapply(1:10, function(i) {
  random.walk.prob(walk_path.RNA.15.list[[i]], mx.coordinate.RNA.15.row.list[[i]])
})) %>% dplyr::mutate(kmer.type = "RNA", kmer.length = 15)

df.prob.initiate.B <- bind_rows(df.prob.DNA.13.initiate.B, df.prob.DNA.15.initiate.B, df.prob.RNA.13.initiate.B, df.prob.RNA.15.initiate.B)

# initiation C
df.prob.DNA.13.initiate.C <- combine.dataframe.list(lapply(1:10, function(i) {
  random.walk.prob(walk_path.DNA.13.list[[i]], mx.coordinate.DNA.13.row.list[[i]])
})) %>% dplyr::mutate(kmer.type = "DNA", kmer.length = 13)

df.prob.DNA.15.initiate.C <- combine.dataframe.list(lapply(1:10, function(i) {
  random.walk.prob(walk_path.DNA.15.list[[i]], mx.coordinate.DNA.15.row.list[[i]])
})) %>% dplyr::mutate(kmer.type = "DNA", kmer.length = 15)

df.prob.RNA.13.initiate.C <- combine.dataframe.list(lapply(1:10, function(i) {
  random.walk.prob(walk_path.RNA.13.list[[i]], mx.coordinate.RNA.13.row.list[[i]])
})) %>% dplyr::mutate(kmer.type = "RNA", kmer.length = 13)

df.prob.RNA.15.initiate.C <- combine.dataframe.list(lapply(1:10, function(i) {
  random.walk.prob(walk_path.RNA.15.list[[i]], mx.coordinate.RNA.15.row.list[[i]])
})) %>% dplyr::mutate(kmer.type = "RNA", kmer.length = 15)

df.prob.initiate.C <- bind_rows(df.prob.DNA.13.initiate.C, df.prob.DNA.15.initiate.C, df.prob.RNA.13.initiate.C, df.prob.RNA.15.initiate.C)

# initiation D
df.prob.DNA.13.initiate.D <- combine.dataframe.list(lapply(1:10, function(i) {
  random.walk.prob(walk_path.DNA.13.list[[i]], mx.coordinate.DNA.13.row.list[[i]])
})) %>% dplyr::mutate(kmer.type = "DNA", kmer.length = 13)

df.prob.DNA.15.initiate.D <- combine.dataframe.list(lapply(1:10, function(i) {
  random.walk.prob(walk_path.DNA.15.list[[i]], mx.coordinate.DNA.15.row.list[[i]])
})) %>% dplyr::mutate(kmer.type = "DNA", kmer.length = 15)

df.prob.RNA.13.initiate.D <- combine.dataframe.list(lapply(1:10, function(i) {
  random.walk.prob(walk_path.RNA.13.list[[i]], mx.coordinate.RNA.13.row.list[[i]])
})) %>% dplyr::mutate(kmer.type = "RNA", kmer.length = 13)

df.prob.RNA.15.initiate.D <- combine.dataframe.list(lapply(1:10, function(i) {
  random.walk.prob(walk_path.RNA.15.list[[i]], mx.coordinate.RNA.15.row.list[[i]])
})) %>% dplyr::mutate(kmer.type = "RNA", kmer.length = 15)

df.prob.initiate.D <- bind_rows(df.prob.DNA.13.initiate.D, df.prob.DNA.15.initiate.D, df.prob.RNA.13.initiate.D, df.prob.RNA.15.initiate.D)

# initiation R
df.prob.DNA.13.initiate.R <- combine.dataframe.list(lapply(1:10, function(i) {
  random.walk.prob(walk_path.DNA.13.list[[i]], mx.coordinate.DNA.13.row.list[[i]])
})) %>% dplyr::mutate(kmer.type = "DNA", kmer.length = 13)

df.prob.DNA.15.initiate.R <- combine.dataframe.list(lapply(1:10, function(i) {
  random.walk.prob(walk_path.DNA.15.list[[i]], mx.coordinate.DNA.15.row.list[[i]])
})) %>% dplyr::mutate(kmer.type = "DNA", kmer.length = 15)

df.prob.RNA.13.initiate.R <- combine.dataframe.list(lapply(1:10, function(i) {
  random.walk.prob(walk_path.RNA.13.list[[i]], mx.coordinate.RNA.13.row.list[[i]])
})) %>% dplyr::mutate(kmer.type = "RNA", kmer.length = 13)

df.prob.RNA.15.initiate.R <- combine.dataframe.list(lapply(1:10, function(i) {
  random.walk.prob(walk_path.RNA.15.list[[i]], mx.coordinate.RNA.15.row.list[[i]])
})) %>% dplyr::mutate(kmer.type = "RNA", kmer.length = 15)

df.prob.initiate.R <- bind_rows(df.prob.DNA.13.initiate.R, df.prob.DNA.15.initiate.R, df.prob.RNA.13.initiate.R, df.prob.RNA.15.initiate.R)
