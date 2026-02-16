Retrieve.unique.kmer.wright <- function(df.dna, df.rna, df.kmer.weight) {
  df.kmer.weight.DNA <- df.kmer.weight %>% dplyr::filter(kmer.type == "DNA")
  df.kmer.weight.RNA <- df.kmer.weight %>% dplyr::filter(kmer.type == "RNA")
  
  df.dna.unq <- dplyr::anti_join(df.dna, df.rna, by = "kmer") %>% dplyr::select("kmer", "RMSE") %>% dplyr::mutate(uniq = "DNA.kmer")
  
  df.dna.unq.count <- merge(df.dna.unq, df.kmer.weight.DNA, by = "kmer")
  
  df.rna.unq <- dplyr::anti_join(df.rna, df.dna, by = "kmer") %>% dplyr::select("kmer", "RMSE") %>% dplyr::mutate(uniq = "RNA.kmer")
  
  df.rna.unq.count <- merge(df.rna.unq, df.kmer.weight.RNA, by = "kmer")
  
  df.out <- dplyr::bind_rows(df.dna.unq.count, df.rna.unq.count)
  
  return(df.out)
}

Unique.kmer.weight.subtypes <- function(df.dna, df.rna, df.kmer.weight) {
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
  
  df.A <- Retrieve.unique.kmer.wright(df.dna.A, df.rna.A, df.kmer.weight) %>% dplyr::mutate(group = "A_enriched")
  df.B <- Retrieve.unique.kmer.wright(df.dna.B, df.rna.B, df.kmer.weight) %>% dplyr::mutate(group = "B_enriched")
  df.C <- Retrieve.unique.kmer.wright(df.dna.C, df.rna.C, df.kmer.weight) %>% dplyr::mutate(group = "C_enriched")
  df.D <- Retrieve.unique.kmer.wright(df.dna.D, df.rna.D, df.kmer.weight) %>% dplyr::mutate(group = "D_enriched")
  df.Rest <- Retrieve.unique.kmer.wright(df.dna.Rest, df.rna.Rest, df.kmer.weight) %>% dplyr::mutate(group = "rest_enriched")
  
  df.fi <- dplyr::bind_rows(df.A, df.B, df.C, df.D, df.Rest) 
  
  return(df.fi)
}
