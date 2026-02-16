mapped.single.genes.dna <- function(df.dna, df.rna, df.map) {
  df.dna <- df.dna %>% dplyr::select(kmer) %>% unique()
  df.rna <- df.rna %>% dplyr::select(kmer) %>% unique()
  
  all.unique.kmer <- dplyr::bind_rows(df.dna, df.rna) %>% dplyr::select(kmer) %>% unique()
  Kmer.dna.uniq <- dplyr::anti_join(df.dna, df.rna, by = "kmer") %>% unique()
  Kmer.dna.uniq.map <- merge(df.map, Kmer.dna.uniq, by = "kmer") 
  
  return(Kmer.dna.uniq.map)
}

mapped.single.genes.rna <- function(df.dna, df.rna, df.map) {
  df.dna <- df.dna %>% dplyr::select(kmer) %>% unique()
  df.rna <- df.rna %>% dplyr::select(kmer) %>% unique()
  
  all.unique.kmer <- dplyr::bind_rows(df.dna, df.rna) %>% dplyr::select(kmer) %>% unique()
  Kmer.rna.uniq <- dplyr::anti_join(df.rna, df.dna, by = "kmer") %>% unique()
  Kmer.rna.uniq.map <- merge(df.map, Kmer.rna.uniq, by = "kmer") 
  
  return(Kmer.rna.uniq.map)
}

plot.uniq.kmer.heatmap <- function(df.dna.uniq, df.rna.uniq) {
  df.dna.uniq_ <- df.dna.uniq %>% dplyr::mutate(order = case_when(group == "A_enriched" ~ "5",
                                                                  group == "B_enriched" ~ "4",
                                                                  group == "C_enriched" ~ "3",
                                                                  group == "D_enriched" ~ "2",
                                                                  group == "rest_enriched" ~ "1"))
  
  df.rna.uniq_ <- df.rna.uniq %>% dplyr::mutate(order = case_when(group == "A_enriched" ~ "5",
                                                                  group == "B_enriched" ~ "4",
                                                                  group == "C_enriched" ~ "3",
                                                                  group == "D_enriched" ~ "2",
                                                                  group == "rest_enriched" ~ "1"))
  
  df.bi <- dplyr::bind_rows(df.dna.uniq_, df.rna.uniq_) %>% dplyr::filter(group != "conserved")
  
  ggplot(df.bi, aes(x = reference_sequence_position, xend = reference_sequence_position + 10, y = order, yend = order, color = group))+geom_segment(linewidth = 2)+theme_bw()+facet_grid(kmer.type ~ .)+xlab("HIV-1 genome")+ylab("HIV-1 subtypes")+xlim(0, 10000)+scale_color_manual(values = c("green4", "orange", "red3", "purple3", "navy"))+theme(axis.title.x=element_text(size=10), axis.text.x=element_text(size=10, colour = "black",angle = 0, vjust = 0), axis.title.y = element_text(size=10),axis.text.y = element_text(size = 10, colour = "black"))
}
