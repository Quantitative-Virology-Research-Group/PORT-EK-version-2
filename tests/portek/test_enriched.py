import pytest
import os
import yaml
import pickle
import pandas as pd
import numpy as np
from unittest.mock import patch, mock_open
from portek import EnrichedKmersPipeline

class TestEnrichedKmersPipelineInit:
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_valid_init(self, mock_read, mock_load, test_project_dir, valid_config_no_header):
        #Setup
        mock_load.return_value = valid_config_no_header
        mock_read.return_value.seq = "ATGCATGC"
        os.makedirs(test_project_dir / "input/indices")
        with open(test_project_dir / "config.yaml", "w") as f:
            yaml.safe_dump(valid_config_no_header, f)
        kmer_set_path1 = test_project_dir / "input/indices/5mer_group1_set.pkl"
        with open(kmer_set_path1, "wb") as f:
            pickle.dump({"AAAAA", "CCCCC"}, f)
        kmer_set_path1 = test_project_dir / "input/indices/5mer_group2_set.pkl"
        with open(kmer_set_path1, "wb") as f:
            pickle.dump({"GGGGG", "TTTTT"}, f)
        sample_list_path1 = test_project_dir / "input/indices/group1_sample_list.pkl"
        with open(sample_list_path1, "wb") as f:
            pickle.dump(["sample1", "sample2"], f)
        sample_list_path2 = test_project_dir / "input/indices/group2_sample_list.pkl"
        with open(sample_list_path2, "wb") as f:
            pickle.dump(["sample1", "sample2"], f)

        #Execute
        enriched_pipeline = EnrichedKmersPipeline(test_project_dir, k=5)

        #Verify
        assert enriched_pipeline.avg_cols == ["group1_avg","group2_avg"]
        assert enriched_pipeline.freq_cols == ["group1_freq","group2_freq"]
        assert enriched_pipeline.c_cols == ["group1_c","group2_c"]
        assert enriched_pipeline.f_cols == ["group1_f","group2_f"]
        assert set(enriched_pipeline.kmer_set) == {"AAAAA", "CCCCC", "GGGGG", "TTTTT"}
        assert sorted(enriched_pipeline.sample_list) == ["group1_sample1", "group1_sample2", "group2_sample1", "group2_sample2"]


    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_init_no_kmers(self, mock_read, mock_load, test_project_dir, valid_config_no_header):
        #Setup
        mock_load.return_value = valid_config_no_header
        mock_read.return_value.seq = "ATGCATGC"
        os.makedirs(test_project_dir / "input/indices")
        with open(test_project_dir / "config.yaml", "w") as f:
            yaml.safe_dump(valid_config_no_header, f)
        kmer_set_path1 = test_project_dir / "input/indices/5mer_group1_set.pkl"
        with open(kmer_set_path1, "wb") as f:
            pickle.dump({"AAAAA", "CCCCC"}, f)
        kmer_set_path1 = test_project_dir / "input/indices/5mer_group2_set.pkl"
        with open(kmer_set_path1, "wb") as f:
            pickle.dump({"GGGGG", "TTTTT"}, f)
        sample_list_path1 = test_project_dir / "input/indices/group1_sample_list.pkl"
        with open(sample_list_path1, "wb") as f:
            pickle.dump(["sample1", "sample2"], f)
        sample_list_path2 = test_project_dir / "input/indices/group2_sample_list.pkl"
        with open(sample_list_path2, "wb") as f:
            pickle.dump(["sample1", "sample2"], f)

        #Execute & Verify
        with pytest.raises(FileNotFoundError, match="Some or all 11-mers are missing from the project directory! Please run PORTEK find_k!"):
            enriched_pipeline = EnrichedKmersPipeline(test_project_dir, k=11)
