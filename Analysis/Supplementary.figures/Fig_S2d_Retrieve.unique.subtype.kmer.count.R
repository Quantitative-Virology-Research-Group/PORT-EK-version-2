Retrieve.unique.subtype.kmer.count <- function(df.dna, df.rna, df.subtype.kmer.count) {
  df.subtype.kmer.count.DNA <- df.subtype.kmer.count %>% dplyr::filter(type == "DNA")
  df.subtype.kmer.count.RNA <- df.subtype.kmer.count %>% dplyr::filter(type == "RNA")
  
  df.dna.unq <- dplyr::anti_join(df.dna, df.rna, by = "kmer") %>% dplyr::select("kmer") %>% dplyr::mutate(uniq = "DNA.kmer")
  
  df.dna.unq.count <- merge(df.dna.unq, df.subtype.kmer.count.DNA, by = "kmer")
  
  df.rna.unq <- dplyr::anti_join(df.rna, df.dna, by = "kmer") %>% dplyr::select("kmer") %>% dplyr::mutate(uniq = "RNA.kmer")
  
  df.rna.unq.count <- merge(df.rna.unq, df.subtype.kmer.count.RNA, by = "kmer")
  
  df.out <- dplyr::bind_rows(df.dna.unq.count, df.rna.unq.count)
  
  return(df.out)
}

Unique.subtype.kmer.count <- function(df.dna, df.rna, df.subtype.kmer.count) {
  df.dna.A <- df.dna %>% dplyr::filter(group == "A_enriched")
  df.dna.B <- df.dna %>% dplyr::filter(group == "B_enriched")
  df.dna.C <- df.dna %>% dplyr::filter(group == "C_enriched")
  df.dna.D <- df.dna %>% dplyr::filter(group == "D_enriched")
  df.dna.Rest <- df.dna %>% dplyr::filter(group == "rest_enriched")
  
  df.rna.A <- df.rna %>% dplyr::filter(group == "A_enriched")
  df.rna.B <- df.rna %>% dplyr::filter(group == "B_enriched")
  df.rna.C <- df.rna %>% dplyr::filter(group == "C_enriched")
  df.rna.D <- df.rna %>% dplyr::filter(group == "D_enriched")
  df.rna.Rest <- df.rna %>% dplyr::filter(group == "rest_enriched")
  
  df.A <- Retrieve.unique.subtype.kmer.count(df.dna.A, df.rna.A, df.subtype.kmer.count) %>% dplyr::mutate(group = "A_enriched")
  df.B <- Retrieve.unique.subtype.kmer.count(df.dna.B, df.rna.B, df.subtype.kmer.count) %>% dplyr::mutate(group = "B_enriched")
  df.C <- Retrieve.unique.subtype.kmer.count(df.dna.C, df.rna.C, df.subtype.kmer.count) %>% dplyr::mutate(group = "C_enriched")
  df.D <- Retrieve.unique.subtype.kmer.count(df.dna.D, df.rna.D, df.subtype.kmer.count) %>% dplyr::mutate(group = "D_enriched")
  df.Rest <- Retrieve.unique.subtype.kmer.count(df.dna.Rest, df.rna.Rest, df.subtype.kmer.count) %>% dplyr::mutate(group = "rest_enriched")
  
  df.fi <- dplyr::bind_rows(df.A, df.B, df.C, df.D, df.Rest) 
  
  return(df.fi)
}
