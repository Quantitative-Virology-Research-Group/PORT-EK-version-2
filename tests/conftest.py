import pytest
import os
import yaml
import numpy as np
import pandas as pd
from unittest.mock import patch, mock_open, MagicMock
from portek.portek_utils import BasePipeline
from portek.portek_map import MappingPipeline
from portek.portek_enriched import EnrichedKmersPipeline


@pytest.fixture
def valid_config():
    return {
        "sample_groups": ["group1", "group2"],
        "input_files": ["group1.fasta", "group2.fasta"],
        "header_format": ["gisaid", "ncbi"],
        "mode": "ovr",
        "goi": "group1",
        "ref_seq": "ref_seq.fasta",
        "ref_genes": [
            {"gene": "gene1", "start": 1, "end": 1},
            {"gene": "gene2", "start": 2, "end": 2},
        ],
    }


@pytest.fixture
def valid_config_no_header():
    return {
        "sample_groups": ["group1", "group2"],
        "input_files": ["group1.fasta", "group2.fasta"],
        "header_format": [],
        "mode": "ovr",
        "goi": "group1",
        "ref_seq": "ref_seq.fasta",
    }


@pytest.fixture
def invalid_config_lenghts():
    return {
        "sample_groups": ["group1"],
        "input_files": ["group1.fasta", "group2.fasta"],
        "header_format": ["format1"],
        "mode": "invalid_mode",
        "goi": "group1",
        "ref_seq": "ref_seq.fasta",
    }


@pytest.fixture
def invalid_config_headers_lengths():
    return {
        "sample_groups": ["group1", "group2"],
        "input_files": ["group1.fasta", "group2.fasta"],
        "header_format": ["format1"],
        "mode": "invalid_mode",
        "goi": "group1",
        "ref_seq": "ref_seq.fasta",
    }


@pytest.fixture
def invalid_config_headers():
    return {
        "sample_groups": ["group1", "group2"],
        "input_files": ["group1.fasta", "group2.fasta"],
        "header_format": ["unknown", "something"],
        "mode": "ovr",
        "goi": "group1",
        "ref_seq": "ref_seq.fasta",
    }


@pytest.fixture
def invalid_config_mode():
    return {
        "sample_groups": ["group1", "group2"],
        "input_files": ["group1.fasta", "group2.fasta"],
        "header_format": ["unknown", "something"],
        "mode": "invalid_mode",
        "goi": "group1",
        "ref_seq": "ref_seq.fasta",
    }


@pytest.fixture
def test_project_dir(tmp_path):
    return tmp_path


@pytest.fixture
@patch("builtins.open", new_callable=mock_open)
@patch("yaml.safe_load")
@patch("Bio.SeqIO.read")
def test_base_pipeline(mock_read, mock_load, mock_open, test_project_dir, valid_config):
    mock_load.return_value = valid_config
    mock_read.return_value.seq = "AAAAATTTTTCCCCCGGGGG"

    return BasePipeline(str(test_project_dir))


# @pytest.fixture
# def mock_config_ava():
#     return {
#         "sample_groups": ["group1", "group2"],
#         "mode": "ava",
#         "goi": "group1",
#         "ref_seq": "reference.fasta",
#         "ref_genes": ["gene1", "gene2"],
#     }


@pytest.fixture
def test_enriched_kmer_stats_csv():
    return pd.DataFrame(
        {
            "group1_avg": [
                0.0,
                1.0,
                1.0,
            ],
            "group2_avg": [
                1.0,
                0.0,
                1.0,
            ],
            "group1-group2_err": [
                -1.0,
                1.0,
                0.0,
            ],
            "group1-group2_p-value": [
                0.001,
                0.001,
                1.0,
            ],
            "RMSE": [
                1.0,
                1.0,
                0.0,
            ],
            "group": [
                "group2_enriched",
                "group1_enriched",
                "conserved",
            ],
            "exclusivity": [
                "exclusive",
                "exclusive",
                "non-exclusive",
            ],
        },
        index=[
            "TTCGA",
            "TACGA",
            "CGAAA",
        ],
        columns=[
            "group1_avg",
            "group2_avg",
            "group1-group2_err",
            "group1-group2_p-value",
            "RMSE",
            "group",
            "exclusivity",
        ],
    )


# @pytest.fixture
# @patch("os.path.isdir")
# @patch("builtins.open", new_callable=mock_open)
# @patch("yaml.safe_load")
# @patch("pandas.read_csv")
# @patch("Bio.SeqIO.read")
# def mock_proper_mapping_pipeline(
#     mock_seqio_read,
#     mock_read_csv,
#     mock_yaml_load,
#     mock_open,
#     mock_isdir,
#     mock_config_ava,
#     mock_enriched_csv,
# ):
#     mock_isdir.return_value = True
#     mock_yaml_load.return_value = mock_config_ava
#     mock_read_csv.return_value = mock_enriched_csv
#     mock_seqio_read.return_value.seq = "ATGC"

#     return MappingPipeline("/fake/dir", 5)


# @pytest.fixture
# def actual_positions():
#     return {
#         "group1": {"AAAAA": np.array([10, 20, 30]), "CCCCC": np.array([15, 25, 35])},
#         "group2": {"GGGGG": np.array([5, 15, 25]), "TTTTT": np.array([10, 20, 30])},
#     }


# @pytest.fixture
# @patch("os.path.isdir")
# @patch("builtins.open", new_callable=mock_open)
# @patch("yaml.safe_load")
# def mock_enriched_pipeline_ava(
#     mock_yaml_load,
#     mock_open,
#     mock_isdir,
#     mock_config_ava
# ):
#     mock_isdir.return_value = True
#     mock_yaml_load.return_value = mock_config_ava

#     return EnrichedKmersPipeline("/fake/dir", 5)
