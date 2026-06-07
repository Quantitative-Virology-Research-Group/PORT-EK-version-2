## Computation of "*K*-mer weight"

This attribute was designated to unveil the composition of enriched *k*-mers. Here, we borrowed the concept of ordinal encoding, which presents baseline encoding, to represent each of the four nucleotides (A, T, C, and G) as a number between 0 and 1. We first utilized the previous setting of ordinal encoding from the study of Wade et al. (2024), whereby A is encoded as 0.25, T as 0.50, C as 0.75, and G as 1.00, and summed the weights for individual enriched k-mers. The “k-mer weight” is defined by the Equation below.

>$$Kmer\ weight = \sum_{n=1}^{n} c(w, n)$$

, where $n$ is the index of each nucleotide over a string of an enriched *k*-mer sequence, and $c(w,n)$ represents the assigned ordinal encoding value $w$ in a specific nucleotide $n$. 
### R code
```
Ordinal.encoding <- function(df) {
  kmer <- df %>% dplyr::select(kmer)
  max_of_rows <- as.numeric(dim(df)[1])
  
  df.fi <- data.frame()
  
  for(i in 1:nrow(df)) {
    if(i > max_of_rows) {
      break
    } else {
      df.tmp <- as.data.frame(strsplit(kmer[i,], "")[[1]])
      names(df.tmp) <- c("kmer.char")
      df.tmp <- df.tmp %>% dplyr::mutate(weight.chr = case_when(kmer.char == "A" ~ 0.25,
                                                                kmer.char == "T" ~ 0.50,
                                                                kmer.char == "C" ~ 0.75,
                                                                kmer.char == "G" ~ 1))
      df.kmer.weight <- as.data.frame(sum(df.tmp$weight.chr))
      names(df.kmer.weight) <- c("kmer.weight")
    }
    df.fi <- rbind(df.fi, df.kmer.weight) 
  }
  df.fii <- dplyr::bind_cols(kmer, df.fi)
  return(df.fii)
}

kmer.weight <- function(df) {
  df.A <- df %>% dplyr::filter(group == "A_enriched")
  df.B <- df %>% dplyr::filter(group == "B_enriched")
  df.C <- df %>% dplyr::filter(group == "C_enriched")
  df.D <- df %>% dplyr::filter(group == "D_enriched")
  df.rest <- df %>% dplyr::filter(group == "rest_enriched")
  
  kmer.weight.A <- Ordinal.encoding(df.A) %>% dplyr::mutate(group = "A")
  kmer.weight.B <- Ordinal.encoding(df.B) %>% dplyr::mutate(group = "B")
  kmer.weight.C <- Ordinal.encoding(df.C) %>% dplyr::mutate(group = "C")
  kmer.weight.D <- Ordinal.encoding(df.D) %>% dplyr::mutate(group = "D")
  kmer.weight.R <- Ordinal.encoding(df.rest) %>% dplyr::mutate(group = "R")
  
  df.pool <- dplyr::bind_rows(kmer.weight.A, kmer.weight.B, kmer.weight.C, kmer.weight.D, kmer.weight.R)
  
  return(df.pool)
}
```
##

<img src="https://github.com/Quantitative-Virology-Research-Group/PORT-EK-version-2/blob/main/images/Kmer.weight.png">

(**a**)Box plots representing the percentage of “*k*-mer weight” across different groups of HIV-1 subtypes. Facets at the x-axis separate enriched DNA and RNA *k*-mers; facets at the y-axis separate enriched k-mers in length of 13 bp or 15 bp. Boxes marked green, orange, red, purple, and dark blue represent HIV-1 subtype A, B, C, D, and rare subtypes, respectively. Significance levels are denoted as follows: ns for no significance, **p* 0.05, ***p* 0.01, ****p* 0.001, *****p* 0.0001. (**b**) Box plots representing the percentage of “*k*-mer weight” across different groups of HIV-1 subtypes. Facets at the x-axis separate *k*-mers enriched from different groups of HIV-1 subtypes; facets at the y-axis separate enriched *k*-mers in a length of 13 bp or 15 bp. Significance levels are denoted as follows: ns for no significance, **p* 0.05.

>We observed a modest and significant enrichment of enriched DNA and RNA *k*-mers in subtype C compared with *k*-mers enriched in other groups of HIV-1 subtypes while comparing “*k*-mer weight” across each group (**a**). No difference in “*k*-mer weight” between *k*-mers in a length of 13 bp and 15 bp, except for a minor enrichment in enriched DNA *k*-mers (13 bp) in subtypes B (enriched DNA *k*-mers: 7.564, enriched RNA *k*-mers: 7.501) and C (enriched DNA *k*-mers: 7.463, enriched RNA *k*-mers: 7.416) compared with enriched RNA *k*-mers (**b**). 
