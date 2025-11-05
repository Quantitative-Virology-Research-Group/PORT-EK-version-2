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

count.map.single.gene.proportion <- function(df) {
  LTR5 <- dim(dplyr::filter(df, gene == "5'LTR"))[1]/dim(df)[1]
  gag <- dim(dplyr::filter(df, gene == "gag"))[1]/dim(df)[1]
  pol <- dim(dplyr::filter(df, gene == "pol"))[1]/dim(df)[1]
  vif <- dim(dplyr::filter(df, gene == "vif"))[1]/dim(df)[1]
  vpr <- dim(dplyr::filter(df, gene == "vpr"))[1]/dim(df)[1]
  tat <- dim(dplyr::filter(df, gene == "tat"))[1]/dim(df)[1]
  rev <- dim(dplyr::filter(df, gene == "rev"))[1]/dim(df)[1]
  vpu <- dim(dplyr::filter(df, gene == "vpu"))[1]/dim(df)[1]
  env <- dim(dplyr::filter(df, gene == "env"))[1]/dim(df)[1]
  nef <- dim(dplyr::filter(df, gene == "nef"))[1]/dim(df)[1]
  LTR3 <- dim(dplyr::filter(df, gene == "3'LTR"))[1]/dim(df)[1]
  
  df.proportion <- data.frame(gene = c("5'LTR", "gag", "pol", "vif", "vpr", "tat", "rev", "vpu", "env", "nef", "3'LTR"), value = c(LTR5, gag, pol, vif, vpr, tat, rev, vpu, env, nef, LTR3))
  
  return(df.proportion)
}

trim.map.single.gene.proportion <- function(df) {
  df.hiv.gene <- data.frame(gene = c("5'LTR", "gag", "pol", "vif", "vpr", "tat", "rev", "vpu", "env", "nef", "3'LTR"))
  df.in <- dplyr::inner_join(df, df.hiv.gene, by = "gene")
  df.out <- dplyr::anti_join(df.hiv.gene, df, by = "gene") %>% dplyr::mutate(value = 0)
  
  df.comp <- dplyr::bind_rows(df.in, df.out)
  
  Richtung.gene <- c("5'LTR", "gag", "pol", "vif", "vpr", "tat", "rev", "vpu", "env", "nef", "3'LTR")
  df.comp$gene <- factor(df.comp$gene, levels = Richtung.gene)
  
  return(df.comp)
}

plot.map.sigle.gene <- function(df.dna.uniq, df.rna.uniq) {
  df.A.dna <- df.dna.uniq %>% dplyr::filter(group == "A_enriched")
  df.A.dna.map <- trim.map.single.gene.proportion(count.map.single.gene.proportion(df.A.dna)) %>% dplyr::rename(A.dna = value)
  
  df.B.dna <- df.dna.uniq %>% dplyr::filter(group == "B_enriched")
  df.B.dna.map <- trim.map.single.gene.proportion(count.map.single.gene.proportion(df.B.dna)) %>% dplyr::rename(B.dna = value)
  
  df.C.dna <- df.dna.uniq %>% dplyr::filter(group == "C_enriched")
  df.C.dna.map <- trim.map.single.gene.proportion(count.map.single.gene.proportion(df.C.dna)) %>% dplyr::rename(C.dna = value)
  
  df.D.dna <- df.dna.uniq %>% dplyr::filter(group == "D_enriched")
  df.D.dna.map <- trim.map.single.gene.proportion(count.map.single.gene.proportion(df.D.dna)) %>% dplyr::rename(D.dna = value)
  
  df.R.dna <- df.dna.uniq %>% dplyr::filter(group == "rest_enriched")
  df.R.dna.map <- trim.map.single.gene.proportion(count.map.single.gene.proportion(df.R.dna)) %>% dplyr::rename(R.dna = value)
  
  df.A.rna <- df.rna.uniq %>% dplyr::filter(group == "A_enriched")
  df.A.rna.map <- trim.map.single.gene.proportion(count.map.single.gene.proportion(df.A.rna)) %>% dplyr::rename(A.rna = value)
  
  df.B.rna <- df.rna.uniq %>% dplyr::filter(group == "B_enriched")
  df.B.rna.map <- trim.map.single.gene.proportion(count.map.single.gene.proportion(df.B.rna)) %>% dplyr::rename(B.rna = value)
  
  df.C.rna <- df.rna.uniq %>% dplyr::filter(group == "C_enriched")
  df.C.rna.map <- trim.map.single.gene.proportion(count.map.single.gene.proportion(df.C.rna)) %>% dplyr::rename(C.rna = value)
  
  df.D.rna <- df.rna.uniq %>% dplyr::filter(group == "D_enriched")
  df.D.rna.map <- trim.map.single.gene.proportion(count.map.single.gene.proportion(df.D.rna)) %>% dplyr::rename(D.rna = value)
  
  df.R.rna <- df.rna.uniq %>% dplyr::filter(group == "rest_enriched")
  df.R.rna.map <- trim.map.single.gene.proportion(count.map.single.gene.proportion(df.R.rna)) %>% dplyr::rename(R.rna = value)
  
  
  df.combine <- dplyr::bind_cols(df.A.dna.map, df.B.dna.map, df.C.dna.map, df.D.dna.map, df.R.dna.map, df.A.rna.map, df.B.rna.map, df.C.rna.map, df.D.rna.map, df.R.rna.map) %>% dplyr::select(gene...1, A.dna, B.dna, C.dna, D.dna, R.dna, A.rna, B.rna, C.rna, D.rna, R.rna) %>% dplyr::rename(gene = gene...1)
  
  rownames(df.combine) <- df.combine$gene
  df.combine$gene <- NULL
  df.combine.mx <- data.matrix(df.combine, rownames.force = T)
  
  col_fun <- colorRamp2(c(0, .1, .2, .3, .4), c("midnightblue", "lightblue", "lightyellow", "yellow", "orangered"))
  
  Heatmap(df.combine.mx, col = col_fun, clustering_distance_columns = "euclidean", clustering_distance_rows = "euclidean", show_row_names = T, show_column_names = T)
}
