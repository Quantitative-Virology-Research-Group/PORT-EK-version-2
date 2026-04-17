import pytest
import pandas as pd
from pandas.testing import assert_frame_equal


@pytest.fixture
def expected_DNA_df():
    return pd.read_csv(
        "tests/big_data_tests/uat_reference/HIV_M_DNA/enriched_15mers_stats.csv",
        index_col="kmer",
    )


@pytest.fixture
def expected_RNA_df():
    return pd.read_csv(
        "tests/big_data_tests/uat_reference/HIV_M_RNA/enriched_15mers_stats.csv",
        index_col="kmer",
    )


@pytest.fixture
def actual_DNA_df():
    return pd.read_csv(
        "projects/HIV_M_DNA/output/enriched_15mers_stats.csv",
        index_col="kmer",
    )


@pytest.fixture
def actual_RNA_df():
    return pd.read_csv(
        "projects/HIV_M_RNA/output/enriched_15mers_stats.csv",
        index_col="kmer",
    )


def test_DNA_df_equivalence(expected_DNA_df, actual_DNA_df):
    assert_frame_equal(expected_DNA_df, actual_DNA_df)


def test_RNA_df_equivalence(expected_RNA_df, actual_RNA_df):
    assert_frame_equal(expected_RNA_df, actual_RNA_df)
