quantile.0_5 <- function(df) {
  df <- df %>% mutate(group_index = case_when( bin < 501 ~ 1,
                                               bin > 500 & bin < 1001 ~ 2,
                                               bin > 1000 & bin < 1501 ~ 3,
                                               bin > 1500 & bin < 2001 ~ 4,
                                               bin > 2000 & bin < 2501 ~ 5,
                                               bin > 2500 & bin < 3001 ~ 6,
                                               bin > 3000 & bin < 3501 ~ 7,
                                               bin > 3500 & bin < 4001 ~ 8,
                                               bin > 4000 & bin < 4501 ~ 9,
                                               bin > 4500 & bin < 5001 ~ 10,
                                               bin > 5000 & bin < 5501 ~ 11,
                                               bin > 5500 & bin < 6001 ~ 12,
                                               bin > 6000 & bin < 6501 ~ 13,
                                               bin > 6500 & bin < 7001 ~ 14,
                                               bin > 7000 & bin < 7501 ~ 15,
                                               bin > 7500 & bin < 8001 ~ 16,
                                               bin > 8000 & bin < 8501 ~ 17,
                                               bin > 8500 & bin < 9001 ~ 18,
                                               bin > 9000 & bin < 9501 ~ 19,
                                               bin > 9500 ~ 20))
  
  df <- df[,-1]
  df.agg <- aggregate(. ~ group_index, df, sum) %>% arrange(group_index)
  
  return(df.agg)
}

quantile.1 <- function(df) {
  df <- df %>% mutate(group_index = case_when( bin < 1001 ~ 1,
                                               bin > 1000 & bin < 2001 ~ 2,
                                               bin > 2000 & bin < 3001 ~ 3,
                                               bin > 3000 & bin < 4001 ~ 4,
                                               bin > 4000 & bin < 5001 ~ 5,
                                               bin > 5000 & bin < 6001 ~ 6,
                                               bin > 6000 & bin < 7001 ~ 7,
                                               bin > 7000 & bin < 8001 ~ 8,
                                               bin > 8000 & bin < 9001 ~ 9,
                                               bin > 9000 ~ 10))
  
  df <- df[,-1]
  df.agg <- aggregate(. ~ group_index, df, sum) %>% arrange(group_index)
  
  return(df.agg)
}

quantile.1_5 <- function(df) {
  df <- df %>% mutate(group_index = case_when( bin < 1501 ~ 1,
                                               bin > 1500 & bin < 3001 ~ 2,
                                               bin > 3000 & bin < 4501 ~ 3,
                                               bin > 4500 & bin < 6001 ~ 4,
                                               bin > 6000 & bin < 7501 ~ 5,
                                               bin > 7500 & bin < 9001 ~ 6,
                                               bin > 9000  ~ 7))
  
  df <- df[,-1]
  df.agg <- aggregate(. ~ group_index, df, sum) %>% arrange(group_index)
  
  return(df.agg)
}

reshape.df <- function(df) {
  df.A <- df %>% dplyr::select(group_index, A_enriched) %>% dplyr::rename(coverage = A_enriched) %>% dplyr::mutate(group = "A")
  df.B <- df %>% dplyr::select(group_index, B_enriched) %>% dplyr::rename(coverage = B_enriched) %>% dplyr::mutate(group = "B")
  df.C <- df %>% dplyr::select(group_index, C_enriched) %>% dplyr::rename(coverage = C_enriched) %>% dplyr::mutate(group = "C")
  df.D <- df %>% dplyr::select(group_index, D_enriched) %>% dplyr::rename(coverage = D_enriched) %>% dplyr::mutate(group = "D")
  df.R <- df %>% dplyr::select(group_index, rest_enriched) %>% dplyr::rename(coverage = rest_enriched) %>% dplyr::mutate(group = "R")
  
  df.pool <- bind_rows(df.A, df.B, df.C, df.D, df.R)
  
  return(df.pool)
}

plot.quantile.coverage <- function(df.dna, df.rna) {
  df.dna.reshape <- reshape.df(df.dna) %>% dplyr::mutate(type = "DNA")
  df.rna.reshape <- reshape.df(df.rna) %>% dplyr::mutate(type = "RNA")
  
  df.combine <- dplyr::bind_rows(df.dna.reshape, df.rna.reshape)
  
  ggplot(df.combine, aes(x = as.factor(group_index), y = coverage, color = group, group = group))+geom_line(linetype = "dashed")+geom_point()+theme_bw()+facet_grid(type ~ .)+xlab("Quantile")+ylab("Coverage count")+scale_color_manual(values = c("green4", "orange", "red3", "purple3", "navy"))+stat_compare_means(comparisons = Vergleichung.group, label = "p.signif")+theme(axis.title.x=element_text(size=12), axis.text.x=element_text(size=12, colour = "black",angle=45,vjust=1, hjust = 1), axis.title.y = element_text(size=12),axis.text.y = element_text(size = 12, colour = "black"))
}
