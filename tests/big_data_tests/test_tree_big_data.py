import pytest
from portek import KmerPhyloTreeConstructor


class TestTreeConstruction:
    def test_benchmark_load_and_format(
        self,
        benchmark,
    ):
        def load_and_format():
            kmer_counts_path = "tests/big_data_tests/uat_reference/HIV_M_DNA/15mer_counts_for_classifier.csv"
            constructor = KmerPhyloTreeConstructor(
                kmer_counts_path, subsample_size=1000
            )
            constructor.format_distance_matrix_for_biopyton()

        benchmark(load_and_format)

    def test_benchmark_construct_tree(
        self,
        benchmark,
    ):
        kmer_counts_path = "tests/big_data_tests/uat_reference/HIV_M_DNA/15mer_counts_for_classifier.csv"
        constructor = KmerPhyloTreeConstructor(kmer_counts_path, subsample_size=1000)
        constructor.format_distance_matrix_for_biopyton()

        benchmark(constructor.construct_tree, method="nj")
