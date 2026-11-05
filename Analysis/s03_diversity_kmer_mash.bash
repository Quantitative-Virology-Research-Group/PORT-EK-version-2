# DNA
$ python scripts/s03_diversity_kmer_mash.py project_rebuttal 15 1500
Mean Mash distance (k=15); diagonal = within-group diversity, off-diagonal = between-group
              A      B      C      D      R  CRF01.AE  CRF02.AG
A         0.080  0.133  0.136  0.132  0.127     0.109     0.104
B         0.133  0.080  0.133  0.106  0.127     0.139     0.134
C         0.136  0.133  0.081  0.133  0.133     0.142     0.138
D         0.132  0.106  0.133  0.076  0.127     0.140     0.133
R         0.127  0.127  0.133  0.127  0.106     0.133     0.122
CRF01.AE  0.109  0.139  0.142  0.140  0.133     0.066     0.117
CRF02.AG  0.104  0.134  0.138  0.133  0.122     0.117     0.071
/home/labadmin/Documents/Arbeitplatz/Projekten/HIV_genome_diverty/HIV-1/2026/Rebuttal/PORT-EK-version-2/scripts/s03_diversity_kmer_mash.py:39: Pandas4Warning: Starting with pandas version 4.0 all arguments of sum will be keyword-only.
  within = pd.Series(np.diag(M), index=G); between = (M.sum(1) - within) / (len(G) - 1)
          within  mean_between
A          0.080         0.123
B          0.080         0.128
C          0.081         0.136
D          0.076         0.129
R          0.106         0.128
CRF01.AE   0.066         0.130
CRF02.AG   0.071         0.124

# RNA
$ python scripts/s03_diversity_kmer_mash.py project_rebuttal_RNA 15 1500
Mean Mash distance (k=15); diagonal = within-group diversity, off-diagonal = between-group
              A      B      C      D      R  CRF01.AE  CRF02.AG
A         0.088  0.134  0.135  0.127  0.129     0.109     0.110
B         0.134  0.082  0.135  0.105  0.132     0.139     0.136
C         0.135  0.135  0.081  0.128  0.134     0.141     0.139
D         0.127  0.105  0.128  0.068  0.127     0.135     0.132
R         0.129  0.132  0.134  0.127  0.111     0.133     0.125
CRF01.AE  0.109  0.139  0.141  0.135  0.133     0.062     0.118
CRF02.AG  0.110  0.136  0.139  0.132  0.125     0.118     0.085
/home/labadmin/Documents/Arbeitplatz/Projekten/HIV_genome_diverty/HIV-1/2026/Rebuttal/PORT-EK-version-2/scripts/s03_diversity_kmer_mash.py:39: Pandas4Warning: Starting with pandas version 4.0 all arguments of sum will be keyword-only.
  within = pd.Series(np.diag(M), index=G); between = (M.sum(1) - within) / (len(G) - 1)
          within  mean_between
A          0.088         0.124
B          0.082         0.130
C          0.081         0.135
D          0.068         0.126
R          0.111         0.130
CRF01.AE   0.062         0.129
CRF02.AG   0.085         0.127
