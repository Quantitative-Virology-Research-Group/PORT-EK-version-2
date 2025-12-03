calculate.percent.subtypes.for.RMSE.log2FC.sig <- function(df) {
  df.sig.dna.dominate <- df %>% dplyr::filter(sig.log2FC == "sig.dna.dominant")
  df.sig.rna.dominate <- df %>% dplyr::filter(sig.log2FC == "sig.rna.dominant")
  
  dna.dominant.A <- dim(dplyr::filter(df.sig.dna.dominate, group.dna == "A_enriched"))[1]/dim(df.sig.dna.dominate)[1]
  dna.dominant.B <- dim(dplyr::filter(df.sig.dna.dominate, group.dna == "B_enriched"))[1]/dim(df.sig.dna.dominate)[1]
  dna.dominant.C <- dim(dplyr::filter(df.sig.dna.dominate, group.dna == "C_enriched"))[1]/dim(df.sig.dna.dominate)[1]
  dna.dominant.D <- dim(dplyr::filter(df.sig.dna.dominate, group.dna == "D_enriched"))[1]/dim(df.sig.dna.dominate)[1]
  dna.dominant.R <- dim(dplyr::filter(df.sig.dna.dominate, group.dna == "rest_enriched"))[1]/dim(df.sig.dna.dominate)[1]
  
  rna.dominant.A <- dim(dplyr::filter(df.sig.rna.dominate, group.rna == "A_enriched"))[1]/dim(df.sig.rna.dominate)[1]
  rna.dominant.B <- dim(dplyr::filter(df.sig.rna.dominate, group.rna == "B_enriched"))[1]/dim(df.sig.rna.dominate)[1]
  rna.dominant.C <- dim(dplyr::filter(df.sig.rna.dominate, group.rna == "C_enriched"))[1]/dim(df.sig.rna.dominate)[1]
  rna.dominant.D <- dim(dplyr::filter(df.sig.rna.dominate, group.rna == "D_enriched"))[1]/dim(df.sig.rna.dominate)[1]
  rna.dominant.R <- dim(dplyr::filter(df.sig.rna.dominate, group.rna == "rest_enriched"))[1]/dim(df.sig.rna.dominate)[1]
  
  df.pool <- data.frame(percentage = c(dna.dominant.A, dna.dominant.B, dna.dominant.C, dna.dominant.D, dna.dominant.R, rna.dominant.A, rna.dominant.B, rna.dominant.C, rna.dominant.D, rna.dominant.R), group = rep(c("A_enrichd", "B_enriched", "C_enriched", "D_enriched", "rest_enriched"), 2), type = c(rep("DNA.kmer", 5), rep("RNA.kmer", 5)))
  
  Richtung.group <- c("A_enrichd", "B_enriched", "C_enriched", "D_enriched", "rest_enriched")
  
  Richtung.type <- c("DNA.kmer", "RNA.kmer")
  
  df.pool$group <- factor(df.pool$group, levels = Richtung.group)
  df.pool$type <- factor(df.pool$type, levels = Richtung.type)
  
  return(df.pool)
}
