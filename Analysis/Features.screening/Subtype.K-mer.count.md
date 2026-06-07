## Computation of "Subtpye *k*-mer count"

We aggregated the count of individual *k*-mers across all isolates present in the same group of the HIV-1 subtype; values were then normalized by the total number of enriched *k*-mers present in the same group of the HIV-1 subtype. The “subtype *k*-mer count” is defined by the Equation below.

>$$Subtype\ kmer\ count = \frac{1}{K} \sum_{i=1}^{i} c(k, i)$$

, where $K$ is the total number of enriched *k*-mers in an indicated group, $i$ is the index over all isolates in the indicated group, and $c(k,i)$ represents the count of enriched *k*-mer $k$ in an isolate $i$. 
### R code
```
kmer.count <- function(df) {
  n.col <- as.numeric(ncol(df)[1])
  n.col <- n.col-1
  
  df.t <- as.data.frame(t(df[,2:n.col]))
  
  df.t <- df.t %>% mutate_if(is.character, as.numeric) %>% dplyr::mutate(kmer.count = rowSums(.)/nrow(df)) %>% dplyr::select(kmer.count) %>% rownames_to_column("kmer") 
  
  return(df.t)
}

kmer.count.normalized <- function(df) {
  A <- dplyr::filter(df, sample_group == "A")
  B <- dplyr::filter(df, sample_group == "B")
  C <- dplyr::filter(df, sample_group == "C")
  D <- dplyr::filter(df, sample_group == "D")
  rest <- dplyr::filter(df, sample_group == "rest")
  
  df.A <- kmer.count(A) %>% dplyr::mutate(subtype = "A")
  df.B <- kmer.count(B) %>% dplyr::mutate(subtype = "B")
  df.C <- kmer.count(C) %>% dplyr::mutate(subtype = "C")
  df.D <- kmer.count(D) %>% dplyr::mutate(subtype = "D")
  df.rest <- kmer.count(rest) %>% dplyr::mutate(subtype = "rest")
  
  df.fi <- dplyr::bind_rows(df.A, df.B, df.C, df.D, df.rest)
  
  return(df.fi)
}
```
##

<img src="https://github.com/Quantitative-Virology-Research-Group/PORT-EK-version-2/blob/main/images/Subtype.kmer.count.png">

(**a**) Box plots representing the percentage of “subtype *k*-mer count" displayed on a logarithmic scale across different groups of HIV-1 subtypes. Facets at the x-axis separate enriched DNA and RNA *k*-mers; facets at the y-axis separate enriched k-mers in length of 13 bp or 15 bp. Boxes marked green, orange, red, purple, and dark blue represent HIV-1 subtype A, B, C, D, and rare subtypes, respectively. Significance levels are denoted as follows: ns for no significance, ***p* 0.01, ****p* 0.001, *****p* 0.0001. (**b**) Box plots representing the percentage of “subtype *k*-mer count" displayed on a logarithmic scale across different groups of HIV-1 subtypes. Facets at the x-axis separate *k*-mers enriched from different groups of HIV-1 subtypes; facets at the y-axis separate enriched *k*-mers in a length of 13 bp or 15 bp. Significance levels are denoted as follows: *****p* 0.0001.

> A more abundant “subtype *k*-mer count” appeared in subtype D, whereas relatively few counts were in subtype B, irrespective of enriched DNA or RNA *k*-mers and their lengths (**a**). A diverse and varied “subtype *k*-mer count” between enriched DNA and RNA k-mers across different groups of HIV-1 subtypes was also observed (**b**).
