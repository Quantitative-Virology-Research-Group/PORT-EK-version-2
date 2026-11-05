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
lumped diag: [0.995 0.968 0.959 0.985 0.931 0.991 0.914]  mean diag = 0.963 (== mean p_within when groups have equal n; any difference in the paper = different run/aggregation)
connected components: 4 sizes: [304   1  44   1]
stationary mass per group: [np.float64(0.166), np.float64(0.145), np.float64(0.146), np.float64(0.172), np.float64(0.108), np.float64(0.169), np.float64(0.093)]
permutation null mean p_within = 0.140 (sd 0.012); observed = 0.963; empirical p = 0.0010
50-step time-averaged home occupancy = 0.903  (uniform 1/7 = 0.143; stationary-mass baseline = 0.143)
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
lumped diag: [0.981 0.948 0.96  0.982 0.852 0.952 0.911]  mean diag = 0.941 (== mean p_within when groups have equal n; any difference in the paper = different run/aggregation)
connected components: 4 sizes: [325   7   9   9]
stationary mass per group: [np.float64(0.182), np.float64(0.103), np.float64(0.242), np.float64(0.194), np.float64(0.084), np.float64(0.097), np.float64(0.097)]
permutation null mean p_within = 0.140 (sd 0.012); observed = 0.941; empirical p = 0.0010
50-step time-averaged home occupancy = 0.817  (uniform 1/7 = 0.143; stationary-mass baseline = 0.143)
t=2: max|Pl^t - exact| = 0.112
t=5: max|Pl^t - exact| = 0.310
t=10: max|Pl^t - exact| = 0.474
