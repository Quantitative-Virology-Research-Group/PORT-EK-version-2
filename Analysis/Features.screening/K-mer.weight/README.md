## Computation of *K*-mer weight

This attribute was designated to unveil the composition of enriched *k*-mers. Here, we borrowed the concept of ordinal encoding, which presents baseline encoding, to represent each of the four nucleotides (A, T, C, and G) as a number between 0 and 1. We first utilized the previous setting of ordinal encoding from the study of Wade et al. (2024), whereby A is encoded as 0.25, T as 0.50, C as 0.75, and G as 1.00, and summed the weights for individual enriched k-mers. The “k-mer weight” is defined by the Equation below.

>$$Kmer\ weight = \sum_{n=1}^{n} c(w, n)$$

, where $n$ is the index of each nucleotide over a string of an enriched *k*-mer sequence, and $c(w,n)$ represents the assigned ordinal encoding value $w$ in a specific nucleotide $n$. 
### R code



