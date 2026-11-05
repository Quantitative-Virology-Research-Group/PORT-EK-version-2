# DNA
$ python scripts/setup_sub.py DNA 80 --seed 1 --data_dir data/DNA_degapped --ref path_to_ref/HIV.HXB2.sequence.fasta 
---
A 199 -> 80
B 8743 -> 80
C 903 -> 80
D 89 -> 80
Rare 81 -> 80
CRF01.AE 1107 -> 80
CRF02.AG 92 -> 80
project written to proj_DNA_80_seed1
---

#$ for s in $(seq 1 100); do   python scripts/setup_sub.py DNA 80 --seed $s --ref path_to_ref/HIV.HXB2.sequence.fasta; done
$ python scripts/s08_replicates.py DNA -k 15 --n_per_group 78 --replicates 100 --workers 4     --data_dir data/DNA_degapped     --full_stats project_rebuttal/output/enriched_15mers_stats.csv     --max_mem 2.0 --out reps_DNA_k15

---
100 replicates, n per group = 78 (after deduplication)
Enriched k-mers per group (mean, SD, 2.5-97.5 percentile):
               mean     std      2.5%     97.5%       min       max
A           3567.8   261.5    3003.7    4046.2    2767.0    4142.0
B           3046.7    97.0    2877.4    3239.2    2847.0    3375.0
C           4418.4   118.4    4222.3    4678.2    4189.0    4734.0
D           3643.5    36.0    3574.8    3714.1    3546.0    3729.0
Rare        2681.8    59.9    2555.2    2797.6    2529.0    2831.0
CRF01.AE    4725.6   130.3    4481.6    5006.8    4436.0    5082.0
CRF02.AG    3465.1    78.2    3318.2    3620.7    3277.0    3679.0
tested    332891.3  4062.0  325678.7  339189.6  316814.0  342649.0

Recurrence of k-mers across replicates (per group: union size; share enriched in >=50% / >=90% of replicates):
          union   ge50   ge90
group                        
A          6775  0.513  0.342
B          6684  0.434  0.304
C          8816  0.480  0.366
CRF01.AE   9191  0.500  0.376
CRF02.AG   5129  0.666  0.540
D          4711  0.776  0.671
Rare       3360  0.822  0.663

Share of replicate markers also enriched for the SAME group in the full-data run (mean over replicates):
 A           0.0
B           0.0
C           0.0
D           0.0
Rare        0.0
CRF01.AE    0.0
CRF02.AG    0.0
dtype: float64
---

# RNA
$ python scripts/setup_sub.py RNA 80 --seed 1 --data_dir data/RNA_degapped --ref project_rebuttal_RNA/input/HIV.HXB2.sequence.fasta
A 868 -> 80
B 2660 -> 80
C 1695 -> 80
D 136 -> 80
Rare 135 -> 80
CRF01.AE 1320 -> 80
CRF02.AG 149 -> 80
project written to proj_RNA_80_seed1

#$ for s in $(seq 1 100); do python scripts/setup_sub.py RNA 80 --seed $s --ref project_rebuttal_RNA/input/HIV.HXB2.sequence.fasta
$ python scripts/s08_replicates.py RNA -k 15 --n_per_group 0 --replicates 100 --workers 4 \
    --data_dir data/RNA_degapped \

---
100 replicates, n per group = 135 (after deduplication)
Enriched k-mers per group (mean, SD, 2.5-97.5 percentile):
               mean     std      2.5%     97.5%       min       max
A           2011.0    95.5    1848.8    2230.5    1841.0    2292.0
B           3420.8    95.3    3250.0    3570.5    3162.0    3707.0
C           4655.2    90.2    4494.0    4830.4    4451.0    4856.0
D           5164.1    42.6    5063.5    5245.7    5058.0    5253.0
Rare        3433.9    25.2    3387.4    3479.0    3369.0    3505.0
CRF01.AE    5392.6    91.5    5232.3    5586.5    5201.0    5680.0
CRF02.AG    2461.7    36.1    2387.2    2534.0    2369.0    2553.0
tested    470260.2  3927.9  462315.6  476470.4  461647.0  478917.0

Recurrence of k-mers across replicates (per group: union size; share enriched in >=50% / >=90% of replicates):
          union   ge50   ge90
group                        
A          4786  0.393  0.259
B          5932  0.562  0.436
C          7371  0.627  0.509
CRF01.AE   8329  0.640  0.517
CRF02.AG   3350  0.726  0.623
D          6050  0.858  0.783
Rare       3928  0.880  0.804

Share of replicate markers also enriched for the SAME group in the full-data run (mean over replicates):
 A           0.0
B           0.0
C           0.0
D           0.0
Rare        0.0
CRF01.AE    0.0
CRF02.AG    0.0
dtype: float64
---
