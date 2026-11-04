# DNA sequences
$ python scripts/s06_cohort_summary.py data/DNA DNA out_dna
             n  unique_seq  countries  year_min  year_max  submissions_proxy  nt_len  gap_frac  amb_frac  exact_dups              top_country
group                                                                                                                                        
A          199         198         20    1985.0    2021.0                 47  8841.0     0.151       0.0           1      RU:46; KE:29; CY:28
B         8743        6801         36    1982.0    2023.0                530  8993.0     0.590       0.0        1942  US:6109; CA:599; AU:510
C          903         871         25    1986.0    2019.0                 83  9006.0     0.408       0.0          32    BW:408; ZA:280; ZM:66
D           89          78         12    1983.0    2017.0                 24  8831.0     0.138       0.0          11        UG:59; KR:6; TD:4
Rare        81          80         21    1983.0    2016.0                 39  8952.0     0.105       0.0           1       BR:20; NG:10; CM:6
CRF01.AE  1107        1040         12    1990.0    2024.0                 79  8999.0     0.371       0.0          67    TH:713; CN:149; JP:99
CRF02.AG    92          92         11    1991.0    2019.0                 25  8830.0     0.152       0.0           0      CM:32; GH:14; FR:12

# RNA sequences
$ python scripts/s06_cohort_summary.py data/RNA RNA out_rna
             n  unique_seq  countries  year_min  year_max  submissions_proxy  nt_len  gap_frac  amb_frac  exact_dups             top_country
group                                                                                                                                       
A          868         847         19    1986.0    2021.0                 83  8834.0     0.226     0.000          21  RW:292; KE:150; RU:139
B         2660        2574         31    1978.0    2020.0                171  8959.0     0.355     0.000          86  US:1766; GB:266; PE:71
C         1695        1652         25    1986.0    2023.0                124  8962.0     0.270     0.000          43  ZA:807; ZM:334; TZ:168
D          136         136         11    1986.0    2020.0                 25  8950.5     0.188     0.001           0      UG:88; KE:31; CD:4
Rare       135         135         18    1987.0    2021.0                 49  8919.0     0.121     0.000           0     ES:38; CM:28; NG:18
CRF01.AE  1320        1291         14    1993.0    2019.0                 91  8736.0     0.278     0.000          29  TH:1029; CN:176; VN:67
CRF02.AG   149         149         19    1999.0    2022.0                 46  8741.0     0.199     0.000           0     GH:43; CM:30; NG:16
