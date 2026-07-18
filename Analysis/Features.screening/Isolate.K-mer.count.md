## Computation of "Isolate *k*-mer count"

We aggregated *k*-mer counts by HIV-1 isolates present in the same group of HIV-1 subtype. The sum was normalized by the total number of isolates present in the same group of the HIV-1 subtype. The “isolate *k*-mer count” is defined by the Equation below.

>$$Isolate\ kmer\ count = \frac{1}{N} \sum_{k=1}^{k} C$$

, where $N$ is the total number of isolates in an indicated group, and $C$ is the summation of all enriched *k*-mer counts (to the $k$-th enriched *k*-mer) in an individual isolate. 
### R code
```
kmer.count.w.isolates <- function(df) {
  
  n.col <- as.numeric(ncol(df)[1])
  
  isolate <- as.data.frame(df[,1])
  names(isolate) <- c("isolate")
  
  group <- as.data.frame(df[,n.col])
  names(group) <- c("group")
  
  df.t <- as.data.frame(df[,2:(n.col-1)])
  
  df.t <- df.t %>% mutate_if(is.character, as.numeric) %>% dplyr::mutate(kmer.count = rowSums(.)/nrow(df)) %>% dplyr::select(kmer.count) 
  
  df.fi <- bind_cols(isolate, group, df.t)
  
  return(df.fi)
}

genome.into.subtype.isolate <- function(df) {
  A <- dplyr::filter(df, sample_group == "A")
  B <- dplyr::filter(df, sample_group == "B")
  C <- dplyr::filter(df, sample_group == "C")
  D <- dplyr::filter(df, sample_group == "D")
  rest <- dplyr::filter(df, sample_group == "rest")
  
  df.A <- kmer.count.w.isolates(A) %>% dplyr::mutate(subtype = "A")
  df.B <- kmer.count.w.isolates(B) %>% dplyr::mutate(subtype = "B")
  df.C <- kmer.count.w.isolates(C) %>% dplyr::mutate(subtype = "C")
  df.D <- kmer.count.w.isolates(D) %>% dplyr::mutate(subtype = "D")
  df.rest <- kmer.count.w.isolates(rest) %>% dplyr::mutate(subtype = "rest")
  
  df.fi <- dplyr::bind_rows(df.A, df.B, df.C, df.D, df.rest)
  
  return(df.fi)
}
```
##

<img src="https://github.com/Quantitative-Virology-Research-Group/PORT-EK-version-2/blob/main/images/Isolate.kmer.count.png">

(**a**) Box plots representing the percentage of “isolate *k*-mer count" across different groups of HIV-1 subtypes. Facets at the x-axis separate enriched DNA and RNA *k*-mers; facets at the y-axis separate enriched k-mers in length of 13 bp or 15 bp. Boxes marked green, orange, red, purple, and dark blue represent HIV-1 subtype A, B, C, D, and rare subtypes, respectively. Significance levels are denoted as follows: ns for no significance, *****p* 0.0001. (**b**) Box plots representing the percentage of “isolate *k*-mer count" across different groups of HIV-1 subtypes. Facets at the x-axis separate *k*-mers enriched from different groups of HIV-1 subtypes; facets at the y-axis separate enriched *k*-mers in a length of 13 bp or 15 bp. Significance levels are denoted as follows: **p* 0.05, ***p* 0.01, ****p* 0.001, *****p* 0.0001.

>While computing “isolate *k*-mer count”, a superior enrichment of enriched DNA k-mers appeared in the group of rare subtypes, followed by subtype D and A (**a**), whereas the majority of enriched RNA *k*-mers accumulated in the HIV-1 subtype D, followed by the rare subtypes and a relatively minor number of *k*-mers were enriched in isolates across the subtypes A, B, and C (**a**). A higher “isolate *k*-mer count” in enriched DNA *k*-mers than in enriched RNA *k*-mers was observed in subtypes A, C, and the rare subtypes (**b**), whereas an opposite pattern was detected in subtypes B and D (**b**).
