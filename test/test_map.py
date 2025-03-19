import pathlib
import pickle
import os
import yaml
from unittest.mock import MagicMock, mock_open, patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from portek.portek_map import MappingPipeline


class TestMappingPipelineInit:
    @patch("builtins.open", new_callable=mock_open)
    @patch("pandas.read_csv")
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_valid_init(
        self,
        mock_read_seq,
        mock_load,
        mock_read_csv,
        mock_open,
        test_project_dir,
        mock_enriched_csv,
        valid_config_no_header,
    ):
        # Setup
        mock_load.return_value = valid_config_no_header
        mock_read_seq.return_value.seq = "ATGCATGC"
        mock_read_csv.return_value = mock_enriched_csv
        # Execute
        mapping_pipeline = MappingPipeline(test_project_dir, 5)

        # Verify
        pd.testing.assert_frame_equal(
            mapping_pipeline.matrices["enriched"], mock_enriched_csv
        )

    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_missing_enriched_kmer_matrix(
        self,
        mock_read_seq,
        mock_load,
        mock_open,
        test_project_dir,
        valid_config_no_header,
    ):
        # Setup
        mock_load.return_value = valid_config_no_header
        mock_read_seq.return_value.seq = "ATGCATGC"

        # Execute & Verify
        with pytest.raises(
            FileNotFoundError,
            match=f"No enriched 5-mers table found in {test_project_dir}output/ ! Please run PORT-EK find_enriched first!",
        ):
            mapping_pipeline = MappingPipeline(test_project_dir, 5)


class TestMappingPipelineCheckIndex:
    @patch("pandas.read_csv")
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_index_exists(
        self,
        mock_read_seq,
        mock_load,
        mock_read_csv,
        test_project_dir,
        mock_enriched_csv,
        valid_config_no_header,
    ):
        # Setup
        mock_load.return_value = valid_config_no_header
        mock_read_seq.return_value.seq = "ATGCATGC"
        mock_read_csv.return_value = mock_enriched_csv
        os.makedirs(test_project_dir / "temp/ref_index")
        with open(test_project_dir / "config.yaml", "w") as f:
            yaml.safe_dump(valid_config_no_header, f)
        with open(test_project_dir / "temp/ref_index/ref_seq_index_5_1.pkl", "w") as f:
            f.write("")
        with open(test_project_dir / "temp/enriched_5mers.pkl", "wb") as f:
            pickle.dump([0,1,2,3], f)

        # Execute
        mapping_pipeline = MappingPipeline(test_project_dir, 5)
        result = mapping_pipeline._check_index(1)
        # Verify
        assert result == True

    @patch("pandas.read_csv")
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_index_not_exists(
        self,
        mock_read_seq,
        mock_load,
        mock_read_csv,
        test_project_dir,
        mock_enriched_csv,
        valid_config_no_header,
    ):
        # Setup
        mock_load.return_value = valid_config_no_header
        mock_read_seq.return_value.seq = "ATGCATGC"
        mock_read_csv.return_value = mock_enriched_csv
        os.makedirs(test_project_dir / "temp")
        with open(test_project_dir / "config.yaml", "w") as f:
            yaml.safe_dump(valid_config_no_header, f)
        with open(test_project_dir / "temp/enriched_5mers.pkl", "wb") as f:
            pickle.dump([0,1,2,3], f)

        # Execute
        mapping_pipeline = MappingPipeline(test_project_dir, 5)
        result = mapping_pipeline._check_index(1)
        # Verify
        assert result == False


class TestMappingPipelineIndexRefSeq:
    @patch("pickle.load")
    @patch("portek.portek_map.MappingPipeline._check_index")
    @patch("pandas.read_csv")
    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_index_exists(
        self,
        mock_read_seq,
        mock_load_yaml,
        mock_open,
        mock_read_csv,
        mock_check_index,
        mock_load_pickle,
        test_project_dir,
        mock_enriched_csv,
        valid_config_no_header,
    ):
        # Setup
        mock_load_yaml.return_value = valid_config_no_header
        mock_read_seq.return_value.seq = "ATGCATGC"
        mock_read_csv.return_value = mock_enriched_csv
        mock_check_index.return_value = True
        mock_load_pickle.return_value = "dupa"

        # Execute
        mapping_pipeline = MappingPipeline(test_project_dir, 5)
        mapping_pipeline.index_ref_seq(1)

        # Verify
        assert mapping_pipeline.kmer_index == "dupa"


    @patch("portek.portek_map.MappingPipeline._check_index")
    @patch("pandas.read_csv")
    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_index_not_exists(
        self,
        mock_read_seq,
        mock_load_yaml,
        mock_open,
        mock_read_csv,
        mock_check_index,
        test_project_dir,
        mock_enriched_csv,
        valid_config_no_header,
    ):
        # Setup
        mock_load_yaml.return_value = valid_config_no_header
        mock_read_seq.return_value.seq = "AAAAAAAAAA"
        mock_read_csv.return_value = mock_enriched_csv
        mock_check_index.return_value = False

        expected_kmer_index = {0: {0: [1, 2, 3, 4, 5, 6]}, 1: {0: [1, 2, 3, 4, 5, 6]}}
        # Execute
        mapping_pipeline = MappingPipeline(test_project_dir, 5)
        mapping_pipeline.index_ref_seq(1)

        # Verify
        assert mapping_pipeline.kmer_index == expected_kmer_index


    @patch("portek.portek_map.MappingPipeline._check_index")
    @patch("pandas.read_csv")
    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_complex_ref_seq(
        self,
        mock_read_seq,
        mock_load_yaml,
        mock_open,
        mock_read_csv,
        mock_check_index,
        test_project_dir,
        mock_enriched_csv,
        valid_config_no_header,
    ):
        # Setup
        mock_load_yaml.return_value = valid_config_no_header
        mock_read_seq.return_value.seq = "AAAACAAAAT"

        #AAAAC: AAAC, AAAA
        # AAACA AACA, AAAA, AAAC
        #  AACAA ACAA, AAAA, AACA
        #   ACAAA CAAA, AAAA, ACAA
        #    CAAAA AAAA, CAAA
        #     AAAAT AAAT, AAAA
        mock_read_csv.return_value = mock_enriched_csv
        mock_check_index.return_value = False
        expected_kmer_index = {0: {1:[1],4:[2],16:[3],64:[4],256:[5],3:[6]}, 1: {0:[1,2,3,4,5,6], 1:[1,2], 4:[2,3], 16:[3,4], 64:[4,5], 3:[6]}}
        # Execute
        mapping_pipeline = MappingPipeline(test_project_dir, 5)
        mapping_pipeline.index_ref_seq(1)

        # Verify
        assert mapping_pipeline.kmer_index == expected_kmer_index