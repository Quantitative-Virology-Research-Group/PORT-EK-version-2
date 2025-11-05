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

plot.map.sigle.gene <- function(df.common, df.map) {
  df <- dplyr::inner_join(df.common, df.map, by = "kmer") %>% dplyr::select(kmer, gene, group.x) %>% dplyr::rename(group = group.x)
  
  df.A <- df %>% dplyr::filter(group == "A_enriched")
  df.A.map <- trim.map.single.gene.proportion(count.map.single.gene.proportion(df.A)) %>% dplyr::rename(A = value)
  
  df.B <- df %>% dplyr::filter(group == "B_enriched")
  df.B.map <- trim.map.single.gene.proportion(count.map.single.gene.proportion(df.B)) %>% dplyr::rename(B = value)
  
  df.C <- df %>% dplyr::filter(group == "C_enriched")
  df.C.map <- trim.map.single.gene.proportion(count.map.single.gene.proportion(df.C)) %>% dplyr::rename(C = value)
  
  df.D <- df %>% dplyr::filter(group == "D_enriched")
  df.D.map <- trim.map.single.gene.proportion(count.map.single.gene.proportion(df.D)) %>% dplyr::rename(D = value)
  
  df.R <- df %>% dplyr::filter(group == "rest_enriched")
  df.R.map <- trim.map.single.gene.proportion(count.map.single.gene.proportion(df.R)) %>% dplyr::rename(R = value)
  
  df.combine <- dplyr::bind_cols(df.A.map, df.B.map, df.C.map, df.D.map, df.R.map) %>% dplyr::select(gene...1, A, B, C, D, R) %>% dplyr::rename(gene = gene...1)
  
  rownames(df.combine) <- df.combine$gene
  df.combine$gene <- NULL
  df.combine.mx <- data.matrix(df.combine, rownames.force = T)
  
  col_fun <- colorRamp2(c(0, .1, .2, .3, .4), c("midnightblue", "lightblue", "lightyellow", "yellow", "orangered"))
  
  Heatmap(df.combine.mx, col = col_fun, clustering_distance_columns = "euclidean", clustering_distance_rows = "euclidean", show_row_names = T, show_column_names = T)
}
