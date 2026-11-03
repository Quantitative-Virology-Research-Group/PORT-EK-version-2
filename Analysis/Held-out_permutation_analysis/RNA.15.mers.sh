$ python scripts/s02_heldout_permutation_fdr.py prepare project_rebuttal_RNA -k 15 --cluster_map RNA_degapped_clusters.tsv.clstr --n_jobs 8

> [!NOTE]
> s02 v2026-11-03 (cluster .clstr + truncated-ID fallback)
A discovery 594 held-out 274
B discovery 2124 held-out 536
C discovery 1302 held-out 393
D discovery 92 held-out 44
R discovery 90 held-out 45
CRF01.AE discovery 1027 held-out 293
CRF02.AG discovery 111 held-out 38
Finding all k-mers of lengths 15 to 15 in 7 files.
Done finding all k-mers!

Imported 2132551 kmers and 5340 samples.
1032600 15-mers passed the entropy filter.
The resulting count matrix would take over 2.0 GB of memory. Removing additional 973938 rare 15-mers.
58662 15-mers remaining.

Getting common 15-mer counts.

Identifying enriched 15-mers.

PORT-EK has found:
6058 CRF01.AE_enriched 15-mers
5737 D_enriched 15-mers
4809 C_enriched 15-mers
3566 B_enriched 15-mers
2431 CRF02.AG_enriched 15-mers
2273 A_enriched 15-mers
2250 R_enriched 15-mers
11 conserved 15-mers
pipeline enriched:
 group
CRF01.AE_enriched    6058
D_enriched           5737
C_enriched           4809
B_enriched           3566
CRF02.AG_enriched    2431
A_enriched           2273
R_enriched           2250
conserved              11
Name: count, dtype: int64

DONE. Discovery project: project_rebuttal_RNA_discall   Held-out sequences: project_rebuttal_RNA_test
Next: python scripts/s02_heldout_permutation_fdr.py analyze project_rebuttal_RNA_discall -k 15
