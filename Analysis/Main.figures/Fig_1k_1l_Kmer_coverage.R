df.reshape <- function(df) {
  df.A <- df %>% dplyr::select(bin, A_enriched) %>% dplyr::rename(count = A_enriched) %>% dplyr::mutate(group = "A")
  df.B <- df %>% dplyr::select(bin, B_enriched) %>% dplyr::rename(count = B_enriched) %>% dplyr::mutate(group = "B")
  df.C <- df %>% dplyr::select(bin, C_enriched) %>% dplyr::rename(count = C_enriched) %>% dplyr::mutate(group = "C")
  df.D <- df %>% dplyr::select(bin, D_enriched) %>% dplyr::rename(count = D_enriched) %>% dplyr::mutate(group = "D")
  df.R <- df %>% dplyr::select(bin, rest_enriched) %>% dplyr::rename(count = rest_enriched) %>% dplyr::mutate(group = "R")
  
  df.combine <- dplyr::bind_rows(df.A, df.B, df.C, df.D, df.R)
}

plot.coverage.plot <- function(df.dna, df.rna) {
  df.dna.coverage <- df.reshape(df.dna) %>% dplyr::mutate(kmer.type = "DNA")
  df.rna.coverage <- df.reshape(df.rna) %>% dplyr::mutate(kmer.type = "RNA")
  
  df.coverage.combine <- dplyr::bind_rows(df.dna.coverage, df.rna.coverage)
  
  ggplot(df.coverage.combine, aes(x = bin, y = count, fill = group))+geom_bar(stat = "identity")+theme_bw()+facet_grid(kmer.type ~ group)+xlab("Genome size (bp)")+ylab("K-mer coverage")+scale_fill_manual(values = c("green4", "orange", "red3", "purple3", "navy"))+theme(axis.title.x=element_text(size=12), axis.text.x=element_text(size=12, colour = "black",angle=45,vjust=1, hjust =1), axis.title.y = element_text(size=12),axis.text.y = element_text(size = 12, colour = "black"))
}

df.reshape.v2 <- function(df.dna, df.rna) {
  df.dna.coverage <- df.reshape(df.dna) %>% dplyr::mutate(kmer.type = "DNA")
  df.rna.coverage <- df.reshape(df.rna) %>% dplyr::mutate(kmer.type = "RNA")
  
  df.coverage.combine <- dplyr::bind_rows(df.dna.coverage, df.rna.coverage)
  
  return(df.coverage.combine)
}
