$ python scripts/s02_heldout_permutation_fdr.py prepare project_rebuttal -k 15 --cluster_map DNA_degapped_clusters.tsv.clstr --n_jobs 8

> [!NOTE]
> s02 v2026-11-03 (cluster .clstr + truncated-ID fallback)
A discovery 151 held-out 48
B discovery 6450 held-out 2293
C discovery 617 held-out 286
D discovery 53 held-out 36
R discovery 60 held-out 21
CRF01.AE discovery 954 held-out 153
CRF02.AG discovery 69 held-out 23
Finding all k-mers of lengths 15 to 15 in 7 files.
Done finding all k-mers!

Imported 2019015 kmers and 8354 samples.
1029158 15-mers passed the entropy filter.
The resulting count matrix would take over 2.0 GB of memory. Removing additional 961053 rare 15-mers.
68105 15-mers remaining.

Getting common 15-mer counts.

Identifying enriched 15-mers.

PORT-EK has found:
5319 CRF01.AE_enriched 15-mers
4880 C_enriched 15-mers
4112 A_enriched 15-mers
3735 D_enriched 15-mers
3538 CRF02.AG_enriched 15-mers
3405 B_enriched 15-mers
2932 R_enriched 15-mers
2 conserved 15-mers
pipeline enriched:
 group
CRF01.AE_enriched    5319
C_enriched           4880
A_enriched           4112
D_enriched           3735
CRF02.AG_enriched    3538
B_enriched           3405
R_enriched           2932
conserved               2
Name: count, dtype: int64

DONE. Discovery project: project_rebuttal_discall   Held-out sequences: project_rebuttal_test
