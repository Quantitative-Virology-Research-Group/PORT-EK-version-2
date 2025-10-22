overlap.kmer.RMSE.ratio <- function(df.dna, df.rna) {
  #kmer.all <- dplyr::bind_rows(df.dna, df.rna) %>% dplyr::select(kmer) %>% unique()
  kmer.dna <- df.dna %>% dplyr::select(kmer, RMSE, group) %>% dplyr::rename(RMSE.dna = RMSE, group.dna = group)
  kmer.rna <- df.rna %>% dplyr::select(kmer, RMSE, group) %>% dplyr::rename(RMSE.rna = RMSE, group.rna = group)
  
  common.kmer <- dplyr::inner_join(kmer.dna, kmer.rna, by = "kmer") %>% dplyr::mutate(label = "common")
  uniq.dna.kmer <- dplyr::anti_join(kmer.dna, kmer.rna, by = "kmer") %>% dplyr::mutate(RMSE.rna = 0, group.rna = "non", label = "uniq.dna")
  uniq.rna.kmer <- dplyr::anti_join(kmer.rna, kmer.dna, by = "kmer") %>% dplyr::mutate(RMSE.dna = 0, group.dna = "non", label = "uniq.rna") %>% dplyr::select(kmer, RMSE.dna, group.dna, RMSE.rna, group.rna, label)
  
  df.all <- dplyr::bind_rows(common.kmer, uniq.dna.kmer, uniq.rna.kmer) %>% dplyr::mutate(ratio.RMSE.dna.vs.rna = log2(RMSE.dna/RMSE.rna), sig = case_when(ratio.RMSE.dna.vs.rna > 1 | ratio.RMSE.dna.vs.rna < -1 ~ "sig", .default = "non-sig")) %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x))
  
  df.all <- df.all %>% dplyr::mutate(sig.log2FC = case_when(ratio.RMSE.dna.vs.rna == 0 & sig == "sig" ~ "non-sig",
                                                            ratio.RMSE.dna.vs.rna > 0 & sig == "sig" ~ "sig.dna.dominant",
                                                            ratio.RMSE.dna.vs.rna < 0 & sig == "sig" ~ "sig.rna.dominant",
                                                            .default = "non-sig"))
  
  return(df.all)
}
