import pytest
from pandas import DataFrame
from numpy import ndarray

from portek import KmerPhyloTreeConstructor


@pytest.fixture
def kmer_phylo_tree_instance() -> KmerPhyloTreeConstructor:
    return KmerPhyloTreeConstructor("tests/test_data/test_counts.csv")


class TestKmerPhyloTreeConstructorInit:
    def test_init(self, kmer_phylo_tree_instance: KmerPhyloTreeConstructor):
        assert kmer_phylo_tree_instance is not None

    def test_kmer_counts_df(self, kmer_phylo_tree_instance: KmerPhyloTreeConstructor):

        expected_shape = (21, 28839)
        assert isinstance(
            kmer_phylo_tree_instance.kmer_counts_df,
            DataFrame,
        )
        assert kmer_phylo_tree_instance.kmer_counts_df.shape == expected_shape
        assert "sample_group" not in kmer_phylo_tree_instance.kmer_counts_df.columns

    def test_subsample_size(self):
        subsample_size = 14
        instance = KmerPhyloTreeConstructor(
            "tests/test_data/test_counts.csv", subsample_size=subsample_size
        )
        n_rest = len(
            [name for name in instance.kmer_counts_df.index if name.startswith("rest_")]
        )
        n_D = len(
            [name for name in instance.kmer_counts_df.index if name.startswith("D_")]
        )
        assert instance.kmer_counts_df.shape[0] == subsample_size
        assert instance.kmer_counts_df.shape[1] == 28839
        assert n_rest + n_D == subsample_size
        assert n_rest == 6
        assert n_D == 8

    def test_subsample_size_balanced(self):
        subsample_size = 14
        instance = KmerPhyloTreeConstructor(
            "tests/test_data/test_counts.csv",
            subsample_size=subsample_size,
            balance_groups=True,
        )
        n_rest = len(
            [name for name in instance.kmer_counts_df.index if name.startswith("rest_")]
        )
        n_D = len(
            [name for name in instance.kmer_counts_df.index if name.startswith("D_")]
        )
        assert instance.kmer_counts_df.shape[0] == subsample_size
        assert instance.kmer_counts_df.shape[1] == 28839
        assert n_rest + n_D == subsample_size
        assert n_rest == 7
        assert n_D == 7


class TestCalculateDistanceMatrix:
    def test_calculate_distance_matrix(
        self, kmer_phylo_tree_instance: KmerPhyloTreeConstructor
    ):
        distance_matrix = kmer_phylo_tree_instance._calculate_distance_matrix()
        expected_shape = (21, 21)
        assert isinstance(distance_matrix, ndarray)
        assert distance_matrix.shape == expected_shape
        assert (distance_matrix >= 0).all()


class TestFormatDistanceMatrixForBiopyton:
    def test_format_distance_matrix_for_biopyton(
        self, kmer_phylo_tree_instance: KmerPhyloTreeConstructor
    ):
        kmer_phylo_tree_instance.format_distance_matrix_for_biopyton()
        biopython_dist = kmer_phylo_tree_instance.distance_matrix
        expected_length = 21
        assert isinstance(biopython_dist, list)
        assert len(biopython_dist) == expected_length
        for i, row in enumerate(biopython_dist):
            assert len(row) == i + 1


class TestConstructTree:
    def test_construct_tree(self, kmer_phylo_tree_instance: KmerPhyloTreeConstructor):
        kmer_phylo_tree_instance.format_distance_matrix_for_biopyton()
        kmer_phylo_tree_instance.construct_tree(method="nj")
        tree = kmer_phylo_tree_instance.tree
        assert tree is not None
        assert tree.is_bifurcating()
