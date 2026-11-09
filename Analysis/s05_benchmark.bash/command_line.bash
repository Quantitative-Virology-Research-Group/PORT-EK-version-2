# DNA
$ bash scripts/s05_benchmark.sh data/DNA_degapped/all_DNA.fasta 15 8 5 project_rebuttal_copy_s05
WARNING: kmc not found in PATH -> skipped
WARNING: kmc not found in PATH -> skipped
WARNING: kmc not found in PATH -> skipped
WARNING: kmc not found in PATH -> skipped
WARNING: kmc not found in PATH -> skipped

Summary of successful runs:
                          n  wall_mean  wall_sd  peak_MB_max
tool      task                                             
jellyfish count          5       1.11     0.18          301
portek    count_index    5     173.75    50.70          405
          find_enriched  5     348.72    72.71         8214
          map            5      65.40     6.98         2172
          tree_n250      5      29.27     2.61         6735
raw table: bench_1793973122.tsv

# RNA
$ bash scripts/s05_benchmark.sh data/RNA_degapped/all_RNA.fasta 15 8 5 project_rebuttal_RNA_copy_s05
WARNING: kmc not found in PATH -> skipped
WARNING: kmc not found in PATH -> skipped
WARNING: kmc not found in PATH -> skipped
WARNING: kmc not found in PATH -> skipped
WARNING: kmc not found in PATH -> skipped

Summary of successful runs:
                          n  wall_mean  wall_sd  peak_MB_max
tool      task                                             
jellyfish count          5       1.41     1.20          302
portek    count_index    5      45.06    14.38          322
          find_enriched  5     194.16    37.16         5116
          map            5      60.69    10.48         2171
          tree_n250      5      18.69     2.73         4046
raw table: bench_1793973128.tsv
