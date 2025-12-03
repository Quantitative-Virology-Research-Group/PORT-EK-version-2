kmer.count.w.isolates <- function(df) {
  
  n.col <- as.numeric(ncol(df)[1])
  
  isolate <- as.data.frame(df[,1])
  names(isolate) <- c("isolate")
  
  group <- as.data.frame(df[,n.col])
  names(group) <- c("group")
  
  df.t <- as.data.frame(df[,2:(n.col-1)])
  
  df.t <- df.t %>% mutate_if(is.character, as.numeric) %>% dplyr::mutate(kmer.count = rowSums(.)/nrow(df)) %>% dplyr::select(kmer.count) 
  
  df.fi <- bind_cols(isolate, group, df.t)
  
  return(df.fi)
}

genome.into.subtype.isolate <- function(df) {
  A <- dplyr::filter(df, sample_group == "A")
  B <- dplyr::filter(df, sample_group == "B")
  C <- dplyr::filter(df, sample_group == "C")
  D <- dplyr::filter(df, sample_group == "D")
  rest <- dplyr::filter(df, sample_group == "rest")
  
  df.A <- kmer.count.w.isolates(A) %>% dplyr::mutate(subtype = "A")
  df.B <- kmer.count.w.isolates(B) %>% dplyr::mutate(subtype = "B")
  df.C <- kmer.count.w.isolates(C) %>% dplyr::mutate(subtype = "C")
  df.D <- kmer.count.w.isolates(D) %>% dplyr::mutate(subtype = "D")
  df.rest <- kmer.count.w.isolates(rest) %>% dplyr::mutate(subtype = "rest")
  
  df.fi <- dplyr::bind_rows(df.A, df.B, df.C, df.D, df.rest)
  
  return(df.fi)
}
