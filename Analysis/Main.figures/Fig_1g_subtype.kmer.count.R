kmer.count <- function(df) {
  n.col <- as.numeric(ncol(df)[1])
  n.col <- n.col-1
  
  df.t <- as.data.frame(t(df[,2:n.col]))
  
  df.t <- df.t %>% mutate_if(is.character, as.numeric) %>% dplyr::mutate(kmer.count = rowSums(.)/nrow(df)) %>% dplyr::select(kmer.count) %>% rownames_to_column("kmer") 
  
  return(df.t)
}

kmer.count.normalized <- function(df) {
  A <- dplyr::filter(df, sample_group == "A")
  B <- dplyr::filter(df, sample_group == "B")
  C <- dplyr::filter(df, sample_group == "C")
  D <- dplyr::filter(df, sample_group == "D")
  rest <- dplyr::filter(df, sample_group == "rest")
  
  df.A <- kmer.count(A) %>% dplyr::mutate(subtype = "A")
  df.B <- kmer.count(B) %>% dplyr::mutate(subtype = "B")
  df.C <- kmer.count(C) %>% dplyr::mutate(subtype = "C")
  df.D <- kmer.count(D) %>% dplyr::mutate(subtype = "D")
  df.rest <- kmer.count(rest) %>% dplyr::mutate(subtype = "rest")
  
  df.fi <- dplyr::bind_rows(df.A, df.B, df.C, df.D, df.rest)
  
  return(df.fi)
}
