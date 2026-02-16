Retrieve.common.subtype.kmer.count <- function(df.subtype.kmer.count) {
  
  df.subtype.kmer.count.dna <- df.subtype.kmer.count %>% dplyr::filter(type == "DNA")
  df.subtype.kmer.count.rna <- df.subtype.kmer.count %>% dplyr::filter(type == "RNA")
  
  df.common.kmer <- dplyr::inner_join(df.subtype.kmer.count.dna, df.subtype.kmer.count.rna, by = "kmer") %>% dplyr::select(kmer, kmer.count.x, kmer.count.y, subtype.x) %>% dplyr::rename(kmer.count.dna = kmer.count.x, kmer.count.rna = kmer.count.y, subtype = subtype.x) 
  
  return(df.common.kmer)
}

Common.subtype.kmer.count <- function(df.subtype.kmer.count) {
  df.subtype.kmer.count.A <- df.subtype.kmer.count %>% dplyr::filter(subtype == "A")
  df.subtype.kmer.count.B <- df.subtype.kmer.count %>% dplyr::filter(subtype == "B")
  df.subtype.kmer.count.C <- df.subtype.kmer.count %>% dplyr::filter(subtype == "Cd")
  df.subtype.kmer.count.D <- df.subtype.kmer.count %>% dplyr::filter(subtype == "D")
  df.subtype.kmer.count.Rest <- df.subtype.kmer.count %>% dplyr::filter(subtype == "rest")
  
  df.A <- Retrieve.common.subtype.kmer.count(df.subtype.kmer.count.A) %>% dplyr::mutate(group = "A_enriched")
  df.B <- Retrieve.common.subtype.kmer.count(df.subtype.kmer.count.B) %>% dplyr::mutate(group = "B_enriched")
  df.C <- Retrieve.common.subtype.kmer.count(df.subtype.kmer.count.C) %>% dplyr::mutate(group = "C_enriched")
  df.D <- Retrieve.common.subtype.kmer.count(df.subtype.kmer.count.D) %>% dplyr::mutate(group = "D_enriched")
  df.Rest <- Retrieve.common.subtype.kmer.count(df.subtype.kmer.count.Rest) %>% dplyr::mutate(group = "rest_enriched")
  
  df.fi <- dplyr::bind_rows(df.A, df.B, df.C, df.D, df.Rest) 
  
  return(df.fi)
}
