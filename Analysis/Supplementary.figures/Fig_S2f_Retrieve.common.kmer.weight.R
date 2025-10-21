Retrieve.common.kmer <- function(df.dna, df.rna, df.kmer.weight) {
  
  df.dna.kmer <- df.dna %>% dplyr::select(kmer) %>% unique()
  df.rna.kmer <- df.rna %>% dplyr::select(kmer) %>% unique()
  all.common.kmer <- dplyr::inner_join(df.dna.kmer, df.rna.kmer, by = "kmer") %>% unique()
  
  
  df.common.kmer <- dplyr::inner_join(df.dna, df.rna, by = "kmer") %>% dplyr::select(kmer, RMSE.x, RMSE.y, group.x) %>% dplyr::rename(RMSE.dna = RMSE.x, RMSE.rna = RMSE.y, group = group.x) 
  
  df.common.kmer.weight <- dplyr::inner_join(df.kmer.weight, all.common.kmer, by = "kmer")
  df.common.kmer.weight.dna <- df.common.kmer.weight %>% dplyr::filter(kmer.type == "DNA") %>% dplyr::select(kmer, kmer.weight) %>% dplyr::rename(kmer.weight.dna = kmer.weight)
  df.common.kmer.weight.rna <- df.common.kmer.weight %>% dplyr::filter(kmer.type == "RNA") %>% dplyr::select(kmer, kmer.weight) %>% dplyr::rename(kmer.weight.rna = kmer.weight)
  df.common.kmer.weight.pool <- merge(df.common.kmer.weight.dna, df.common.kmer.weight.rna, by = "kmer")
  
  df.common.kmer.weight.RMSE <- merge(df.common.kmer, df.common.kmer.weight.pool, by = "kmer")
  
  return(df.common.kmer.weight.RMSE)
}
