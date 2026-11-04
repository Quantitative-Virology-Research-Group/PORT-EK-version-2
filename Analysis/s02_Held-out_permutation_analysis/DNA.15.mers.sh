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

$ python scripts/s02_heldout_permutation_fdr.py analyze project_rebuttal_discall -k 15 --n_perm 100 --n_boot 100
s02 v2026-11-03 (cluster .clstr + truncated-ID fallback)
tested k-mers (after entropy/min_freq filter) = 68105; samples = 8354

Enriched k-mers per group: raw vs Benjamini-Yekutieli
           raw_p<0.01  BY-FDR<0.01
A               4112         3000
B               3405         2553
C               4880         4184
CRF01.AE        5319         4538
CRF02.AG        3538         2298
D               3735         2551
R               2932         1339 
retained after FDR: 73.3%

Permutation null (whole selection, 100 shuffles, k-mer-subsample scaled): mean=5 max=27; observed=28199

Bootstrap stability (n=100) of 2000 random markers: median=1.00; fraction >=0.8: 0.80

Held-out assignment (rule: group with highest fraction of its own markers present)
           A     B    C   D   R  CRF01.AE  CRF02.AG
A         48     0    0   0   0         0         0
B          0  2293    0   0   0         0         0
C          0     0  286   0   0         0         0
D          0     0    0  36   0         0         0
R          0     1    2   0  17         0         1
CRF01.AE   0     0    0   0   0       153         0
CRF02.AG   0     0    0   0   0         0        23
              precision    recall  f1-score   support

           A      1.000     1.000     1.000        48
           B      1.000     1.000     1.000      2293
           C      0.993     1.000     0.997       286
           D      1.000     1.000     1.000        36
           R      1.000     0.810     0.895        21
    CRF01.AE      1.000     1.000     1.000       153
    CRF02.AG      0.958     1.000     0.979        23

    accuracy                          0.999      2860
   macro avg      0.993     0.973     0.981      2860
weighted avg      0.999     0.999     0.999      2860

macro-F1 = 0.981
