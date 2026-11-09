# DNA
$ python scripts/s01_isolate_metric_and_graph.py project_rebuttal 15 50 1000
== Isolate metric by group (median) ==
          N_group  orig_sum_over_N  raw_sum  frac_own_group_markers  frac_all_markers
group                                                                                
A           199.0            9.548   1900.0                   0.239             0.070
B          8743.0            0.248   2168.0                   0.440             0.081
C           903.0            2.963   2676.0                   0.452             0.099
CRF01.AE   1107.0            2.631   2913.0                   0.488             0.108
CRF02.AG     92.0           22.620   2081.0                   0.469             0.077
D            89.0           26.247   2336.0                   0.329             0.086
R            81.0           22.432   1817.0                   0.272             0.067

== Random-walk summary ==
markers=350  nodes=350  mean p_within=0.963  p_across=0.037
group order: ['C', 'B', 'D', 'CRF01.AE', 'A', 'CRF02.AG', 'R']
lumped diag (same order): {'C': np.float64(0.995), 'B': np.float64(0.968), 'D': np.float64(0.959), 'CRF01.AE': np.float64(0.985), 'A': np.float64(0.931), 'CRF02.AG': np.float64(0.991), 'R': np.float64(0.914)}  mean diag = 0.963 (== mean p_within when groups have equal n; any difference in the paper = different run/aggregation)
connected components: 4 sizes: [304   1  44   1]
isolated nodes (no edge weight, e.g. isolates carrying none of the markers): 2; mean p_within excluding them = 0.969
  component 1: size 1, groups: {'D': np.int64(1)}
  component 2: size 44, groups: {'A': np.int64(44)}
  component 3: size 1, groups: {'R': np.int64(1)}
stationary (degree-proportional) mass per group: {'C': 0.166, 'B': 0.145, 'D': 0.146, 'CRF01.AE': 0.172, 'A': 0.108, 'CRF02.AG': 0.169, 'R': 0.093} (meaningful only within connected components)
permutation null mean p_within = 0.140 (sd 0.012); observed = 0.963; empirical p = 0.0010
50-step time-averaged home occupancy = 0.903  (uniform 1/7 = 0.143; the mean stationary mass is always 1/7 and is not a useful baseline)
permutation null of 50-step occupancy (200 label shuffles on the same graph): mean = 0.173 (sd 0.004); observed = 0.903; empirical p = 0.0050
t=2: max|Pl^t - exact| = 0.070
t=5: max|Pl^t - exact| = 0.215
t=10: max|Pl^t - exact| = 0.384

# RNA
$ python scripts/s01_isolate_metric_and_graph.py project_rebuttal_RNA 15 50 1000
== Isolate metric by group (median) ==
          N_group  orig_sum_over_N  raw_sum  frac_own_group_markers  frac_all_markers
group                                                                                
A           868.0            2.513   2181.0                   0.393             0.082
B          2660.0            0.962   2560.0                   0.411             0.096
C          1695.0            1.773   3006.0                   0.446             0.113
CRF01.AE   1320.0            2.739   3616.0                   0.536             0.136
CRF02.AG    149.0           13.933   2076.0                   0.430             0.078
D           136.0           22.765   3096.0                   0.479             0.116
R           135.0           16.215   2189.0                   0.227             0.082

== Random-walk summary ==
markers=350  nodes=350  mean p_within=0.941  p_across=0.059
group order: ['C', 'B', 'D', 'CRF01.AE', 'A', 'CRF02.AG', 'R']
lumped diag (same order): {'C': np.float64(0.981), 'B': np.float64(0.948), 'D': np.float64(0.96), 'CRF01.AE': np.float64(0.982), 'A': np.float64(0.852), 'CRF02.AG': np.float64(0.952), 'R': np.float64(0.911)}  mean diag = 0.941 (== mean p_within when groups have equal n; any difference in the paper = different run/aggregation)
connected components: 4 sizes: [325   7   9   9]
isolated nodes (no edge weight, e.g. isolates carrying none of the markers): 0; mean p_within excluding them = 0.941
  component 1: size 7, groups: {'A': np.int64(7)}
  component 2: size 9, groups: {'A': np.int64(9)}
  component 3: size 9, groups: {'A': np.int64(9)}
stationary (degree-proportional) mass per group: {'C': 0.182, 'B': 0.103, 'D': 0.242, 'CRF01.AE': 0.194, 'A': 0.084, 'CRF02.AG': 0.097, 'R': 0.097} (meaningful only within connected components)
permutation null mean p_within = 0.140 (sd 0.012); observed = 0.941; empirical p = 0.0010
50-step time-averaged home occupancy = 0.817  (uniform 1/7 = 0.143; the mean stationary mass is always 1/7 and is not a useful baseline)
permutation null of 50-step occupancy (200 label shuffles on the same graph): mean = 0.174 (sd 0.004); observed = 0.817; empirical p = 0.0050
t=2: max|Pl^t - exact| = 0.112
t=5: max|Pl^t - exact| = 0.310
t=10: max|Pl^t - exact| = 0.474
