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

$ python scripts/s02_heldout_permutation_fdr.py analyze project_rebuttal_RNA_discall -k 15 --n_perm 100 --n_boot 100
s02 v2026-11-03 (cluster .clstr + truncated-ID fallback)
tested k-mers (after entropy/min_freq filter) = 58662; samples = 5340

Enriched k-mers per group: raw vs Benjamini-Yekutieli
           raw_p<0.01  BY-FDR<0.01
A               2273         1807
B               3566         3286
C               4809         4546
CRF01.AE        6058         5563
CRF02.AG        2431         1982
D               5737         5032
R               2250         1903 
retained after FDR: 88.9%

Permutation null (whole selection, 100 shuffles, k-mer-subsample scaled): mean=2 max=35; observed=27193

Bootstrap stability (n=100) of 2000 random markers: median=1.00; fraction >=0.8: 0.85

Held-out assignment (rule: group with highest fraction of its own markers present)
            A    B    C   D   R  CRF01.AE  CRF02.AG
A         274    0    0   0   0         0         0
B           0  536    0   0   0         0         0
C          10    0  383   0   0         0         0
D           0    0    0  44   0         0         0
R           0    0    0   1  44         0         0
CRF01.AE    0    0    0   0   0       293         0
CRF02.AG    0    0    0   0   0         0        38
              precision    recall  f1-score   support

           A      0.965     1.000     0.982       274
           B      1.000     1.000     1.000       536
           C      1.000     0.975     0.987       393
           D      0.978     1.000     0.989        44
           R      1.000     0.978     0.989        45
    CRF01.AE      1.000     1.000     1.000       293
    CRF02.AG      1.000     1.000     1.000        38

    accuracy                          0.993      1623
   macro avg      0.992     0.993     0.992      1623
weighted avg      0.993     0.993     0.993      1623

macro-F1 = 0.992
