Count.unique.common.DNAvsRNA.Kmers <- function(df.dna, df.rna) {
  df.dna <- df.dna %>% dplyr::select(kmer) %>% unique()
  df.rna <- df.rna %>% dplyr::select(kmer) %>% unique()
  
  all.unique.kmer <- dplyr::bind_rows(df.dna, df.rna) %>% dplyr::select(kmer) %>% unique()
  all.unique.kmer.count <- dim(all.unique.kmer)[1]
  
  Kmer.common <- dim(dplyr::inner_join(df.dna, df.rna, by = "kmer"))[1]/all.unique.kmer.count
  Kmer.dna.uniq <- dim(dplyr::anti_join(df.dna, df.rna, by = "kmer"))[1]/all.unique.kmer.count
  Kmer.rna.uniq <- dim(dplyr::anti_join(df.rna, df.dna, by = "kmer"))[1]/all.unique.kmer.count
  
  percentage <- c(Kmer.common, Kmer.dna.uniq, Kmer.rna.uniq)
  label <- c("Common", "DNA.uniq", "RNA.uniq")
  
  df <- data.frame(percentage, label)
}

Separate.unique.common.DNAvsRNA.Kmers.per.group <- function(df.dna, df.rna) {
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
  
  df.A <- Count.unique.common.DNAvsRNA.Kmers(df.dna.A, df.rna.A) %>% dplyr::mutate(group = "A.enriched")
  df.B <- Count.unique.common.DNAvsRNA.Kmers(df.dna.B, df.rna.B) %>% dplyr::mutate(group = "B.enriched")
  df.C <- Count.unique.common.DNAvsRNA.Kmers(df.dna.C, df.rna.C) %>% dplyr::mutate(group = "C.enriched")
  df.D <- Count.unique.common.DNAvsRNA.Kmers(df.dna.D, df.rna.D) %>% dplyr::mutate(group = "D.enriched")
  df.Rest <- Count.unique.common.DNAvsRNA.Kmers(df.dna.Rest, df.rna.Rest) %>% dplyr::mutate(group = "Rest.enriched")
  
  df.fi <- dplyr::bind_rows(df.A, df.B, df.C, df.D, df.Rest)
  
  return(df.fi)
}
