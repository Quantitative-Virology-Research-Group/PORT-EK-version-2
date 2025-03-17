import pytest
import pickle
import os
import pandas as pd
from scipy import stats
from unittest.mock import patch, mock_open, MagicMock
from portek import calc_kmer_pvalue
from portek.portek_utils import BasePipeline


class TestBasePiplelineInit:

    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_load_and_check_config_valid(self,mock_read, mock_load, mock_open,valid_config, test_project_dir):
        #Setup
        mock_load.return_value = valid_config
        mock_read.return_value = MagicMock(seq="ATGCATGC")

        #Execute
        pipeline = BasePipeline(test_project_dir)

        #Verify
        assert pipeline.sample_groups == valid_config["sample_groups"]
        assert pipeline.input_files == valid_config["input_files"]
        assert pipeline.header_format == valid_config["header_format"]
        assert pipeline.mode == valid_config["mode"]
        assert pipeline.goi == valid_config["goi"]
        assert pipeline.ref_seq_name == "ref_seq"
        assert pipeline.ref_seq == "ATGCATGC"

    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_load_and_check_config_valid_no_header(self,
        mock_read, mock_load, mock_open,valid_config_no_header, test_project_dir
    ):
        #Setup
        mock_load.return_value = valid_config_no_header
        mock_read.return_value = MagicMock(seq="ATGCATGC")

        #Execute
        pipeline = BasePipeline(test_project_dir)
    
        #Verify
        assert (
            pipeline.sample_groups
            == valid_config_no_header["sample_groups"]
        )
        assert (
            pipeline.input_files == valid_config_no_header["input_files"]
        )
        assert pipeline.header_format == ["plain","plain"]
        assert pipeline.mode == valid_config_no_header["mode"]
        assert pipeline.goi == valid_config_no_header["goi"]
        assert pipeline.ref_seq_name == "ref_seq"
        assert pipeline.ref_seq == "ATGCATGC"

    def test_load_and_check_config_missing_config(self, test_project_dir):
        with pytest.raises(FileNotFoundError) as context:
            pipeline = BasePipeline(test_project_dir)
            assert "No config.yaml file found in directory" in str(context.value)

    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_load_and_check_config_invalid_config_lengths(self,mock_read, mock_load, mock_open,invalid_config_lenghts, test_project_dir):
        #Setup
        mock_load.return_value = invalid_config_lenghts
        mock_read.return_value = MagicMock(seq="ATGCATGC")

        #Execute & Verify
        with pytest.raises(ValueError) as context:
            pipeline = BasePipeline(test_project_dir)
            assert "Number of sample groups and input fastas must match!" in str(context.value)

    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_load_and_check_config_invalid_config_headers_lengths(self,mock_read, mock_load, mock_open,invalid_config_headers_lengths, test_project_dir):
        #Setup
        mock_load.return_value = invalid_config_headers_lengths
        mock_read.return_value = MagicMock(seq="ATGCATGC")

        #Execute & Verify
        with pytest.raises(ValueError) as context:
            pipeline = BasePipeline(test_project_dir)
            assert "Number of header formats must be 0 or match number of sample groups!" in str(context.value)

    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_load_and_check_config_invalid_config_headers(self,mock_read, mock_load, mock_open,invalid_config_headers, test_project_dir):
        #Setup
        mock_load.return_value = invalid_config_headers
        mock_read.return_value = MagicMock(seq="ATGCATGC")

        #Execute & Verify
        with pytest.raises(ValueError) as context:
            pipeline = BasePipeline(test_project_dir)
            assert "Header format must be 'gisaid', 'ncbi' or 'plain'!" in str(context.value)

    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_load_and_check_config_invalid_config_mode(self,mock_read, mock_load, mock_open,invalid_config_mode, test_project_dir):
        #Setup
        mock_load.return_value = invalid_config_mode
        mock_read.return_value = MagicMock(seq="ATGCATGC")

        #Execute & Verify
        with pytest.raises(ValueError) as context:
            pipeline = BasePipeline(test_project_dir)
            assert "Unrecognized analysis mode, should by ava or ovr. Check your config file!" in str(context.value)

    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    def test_load_and_check_config_missing_reference(self, mock_load, mock_open,valid_config_no_header, test_project_dir):
        #Setup
        mock_load.return_value = valid_config_no_header

        #Execute & Verify
        with pytest.raises(ValueError) as context:
            pipeline = BasePipeline(test_project_dir)
            assert "Missing reference sequence file or the file has incorrect format!" in str(context.value)

    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    def test_load_and_check_config_missing_reference(self, mock_load, mock_open,valid_config_no_header, test_project_dir):
        #Setup
        mock_load.return_value = valid_config_no_header
        with open(test_project_dir / "input/ref_seq.fasta", "w") as f:
            f.write("")

        #Execute & Verify
        with pytest.raises(ValueError) as context:
            pipeline = BasePipeline(test_project_dir)
            assert "Missing reference sequence file or the file has incorrect format!" in str(context.value)           


class TestBasePipelineLoadKmerSet:
    def test_load_kmer_set_success(self, test_base_pipeline, test_project_dir):
        os.makedirs(test_project_dir / "input/indices")
        kmer_set_path = test_project_dir / "input/indices/5mer_group1_set.pkl"
        with open(kmer_set_path, "wb") as f:
            pickle.dump({"AAAAA", "CCCCC"}, f)
        kmer_set_path = test_project_dir / "input/indices/5mer_group2_set.pkl"
        with open(kmer_set_path, "wb") as f:
            pickle.dump({"GGGGG", "TTTTT"}, f)
        test_base_pipeline._load_kmer_set()
        assert set(test_base_pipeline.kmer_set) == {"AAAAA", "CCCCC", "GGGGG", "TTTTT"}

def test_load_kmer_set_missing_files(test_base_pipeline):
    with pytest.raises(FileNotFoundError, match="Some or all k-mers are missing from the project directory! Please run PORTEK find_k!"):
        test_base_pipeline._load_kmer_set()

def test_load_kmer_set_partial_files(test_base_pipeline, test_project_dir):
    os.makedirs(test_project_dir / "input/indices")
    kmer_set_path = test_project_dir / "input/indices/5mer_group1_set.pkl"
    with open(kmer_set_path, "wb") as f:
        pickle.dump({"AAAAA", "CCCCC"}, f)
    with pytest.raises(FileNotFoundError, match="Some or all k-mers are missing from the project directory! Please run PORTEK find_k!"):
        test_base_pipeline._load_kmer_set()

class TestBasePipelineLoadSampleList:
    def test_load_sample_list_success(self, test_base_pipeline, test_project_dir):
        os.makedirs(test_project_dir / "input/indices")
        sample_list_path = test_project_dir / "input/indices/group1_sample_list.pkl"
        with open(sample_list_path, "wb") as f:
            pickle.dump(["sample1", "sample2"], f)
        sample_list_path = test_project_dir / "input/indices/group2_sample_list.pkl"
        with open(sample_list_path, "wb") as f:
            pickle.dump(["sample1", "sample2"], f)
        test_base_pipeline._load_sample_list()
        assert sorted(test_base_pipeline.sample_list) == ["group1_sample1", "group1_sample2", "group2_sample1", "group2_sample2"]

def test_load_sample_list_missing_files(test_base_pipeline):
    with pytest.raises(FileNotFoundError, match="Some or all samples are missing from the project directory! Please run PORTEK find_k!"):
        test_base_pipeline._load_sample_list()

def test_load_sample_list_partial_files(test_base_pipeline, test_project_dir):
    os.makedirs(test_project_dir / "input/indices")
    sample_list_path = test_project_dir / "input/indices/group1_sample_list.pkl"
    with open(sample_list_path, "wb") as f:
        pickle.dump(["sample1", "sample2"], f)
    with pytest.raises(FileNotFoundError, match="Some or all samples are missing from the project directory! Please run PORTEK find_k!"):
        test_base_pipeline._load_sample_list()

class TestCalcKmerPvalue:
    def test_calc_kmer_pvalue_success(self):
        # Setup
        kmer = "AAAAA"
        first_group = ["sample1", "sample2"]
        sec_group = ["sample3", "sample4"]
        data = {
            "sample1": [1, 2, 3],
            "sample2": [2, 3, 4],
            "sample3": [5, 6, 7],
            "sample4": [6, 7, 8],
        }
        matrix = pd.DataFrame(data, index=["AAAAA", "AAAAC", "AAAAG"])

        # Execute
        pvalue = calc_kmer_pvalue(kmer, first_group, sec_group, matrix)

        # Verify
        expected_pvalue = stats.mannwhitneyu(
            matrix.loc[kmer, first_group], matrix.loc[kmer, sec_group]
        ).pvalue
        assert pvalue == expected_pvalue
