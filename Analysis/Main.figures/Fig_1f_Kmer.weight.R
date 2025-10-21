Ordinal.encoding <- function(df) {
  kmer <- df %>% dplyr::select(kmer)
  max_of_rows <- as.numeric(dim(df)[1])
  
  df.fi <- data.frame()
  
  for(i in 1:nrow(df)) {
    if(i > max_of_rows) {
      break
    } else {
      df.tmp <- as.data.frame(strsplit(kmer[i,], "")[[1]])
      names(df.tmp) <- c("kmer.char")
      df.tmp <- df.tmp %>% dplyr::mutate(weight.chr = case_when(kmer.char == "A" ~ 0.25,
                                                                kmer.char == "T" ~ 0.50,
                                                                kmer.char == "C" ~ 0.75,
                                                                kmer.char == "G" ~ 1))
      df.kmer.weight <- as.data.frame(sum(df.tmp$weight.chr))
      names(df.kmer.weight) <- c("kmer.weight")
    }
    df.fi <- rbind(df.fi, df.kmer.weight) 
  }
  df.fii <- dplyr::bind_cols(kmer, df.fi)
  return(df.fii)
}

kmer.weight <- function(df) {
  df.A <- df %>% dplyr::filter(group == "A_enriched")
  df.B <- df %>% dplyr::filter(group == "B_enriched")
  df.C <- df %>% dplyr::filter(group == "C_enriched")
  df.D <- df %>% dplyr::filter(group == "D_enriched")
  df.rest <- df %>% dplyr::filter(group == "rest_enriched")
  
  kmer.weight.A <- Ordinal.encoding(df.A) %>% dplyr::mutate(group = "A")
  kmer.weight.B <- Ordinal.encoding(df.B) %>% dplyr::mutate(group = "B")
  kmer.weight.C <- Ordinal.encoding(df.C) %>% dplyr::mutate(group = "C")
  kmer.weight.D <- Ordinal.encoding(df.D) %>% dplyr::mutate(group = "D")
  kmer.weight.R <- Ordinal.encoding(df.rest) %>% dplyr::mutate(group = "R")
  
  df.pool <- dplyr::bind_rows(kmer.weight.A, kmer.weight.B, kmer.weight.C, kmer.weight.D, kmer.weight.R)
  
  return(df.pool)
}
