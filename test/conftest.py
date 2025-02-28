import pytest
import numpy as np
import pandas as pd
from unittest.mock import patch, mock_open, MagicMock
from portek.portek_map import MappingPipeline

@pytest.fixture
def mock_config():
    return {
        "sample_groups": ["group1", "group2"],
        "mode": "ovr",
        "goi": "group1",
        "ref_seq": "reference.fasta",
        "ref_genes": ["gene1", "gene2"]
    }

@pytest.fixture
def mock_enriched_csv():
    return pd.DataFrame({
        "group": ["group1", "group2"],
        "value": [1, 2]
    })

@pytest.fixture
@patch("os.path.isdir")
@patch("builtins.open", new_callable=mock_open)
@patch("yaml.safe_load")
@patch("pandas.read_csv")
@patch("Bio.SeqIO.read")
def mock_proper_mapping_pipeline(
        mock_seqio_read,
        mock_read_csv,
        mock_yaml_load,
        mock_open,
        mock_isdir,
        mock_config,
        mock_enriched_csv,
    ):
        mock_isdir.return_value = True
        mock_yaml_load.return_value = mock_config
        mock_read_csv.return_value = mock_enriched_csv
        mock_seqio_read.return_value.seq = "ATGC"

        return MappingPipeline("/fake/dir", 5)

@pytest.fixture
def actual_positions():
    return {
        "group1": {
            "AAAAA": np.array([10, 20, 30]),
            "CCCCC": np.array([15, 25, 35])
        },
        "group2": {
            "GGGGG": np.array([5, 15, 25]),
            "TTTTT": np.array([10, 20, 30])
        }
    }