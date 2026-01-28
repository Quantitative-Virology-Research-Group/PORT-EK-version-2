import pytest
import pickle
import os
from unittest.mock import patch, mock_open, MagicMock
import portek
from portek.portek_utils import BasePipeline


class TestBasePiplelineInit:

    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_load_and_check_config_valid(
        self, mock_read, mock_load, mock_open, valid_config, test_project_dir
    ):
        # Setup
        mock_load.return_value = valid_config
        mock_read.return_value = MagicMock(seq="ATGCATGC")

        # Execute
        pipeline = BasePipeline(test_project_dir)

        # Verify
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
    def test_load_and_check_config_valid_no_header(
        self, mock_read, mock_load, mock_open, valid_config_no_header, test_project_dir
    ):
        # Setup
        mock_load.return_value = valid_config_no_header
        mock_read.return_value = MagicMock(seq="ATGCATGC")

        # Execute
        pipeline = BasePipeline(test_project_dir)

        # Verify
        assert pipeline.sample_groups == valid_config_no_header["sample_groups"]
        assert pipeline.input_files == valid_config_no_header["input_files"]
        assert pipeline.header_format == ["plain", "plain"]
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
    def test_load_and_check_config_invalid_config_lengths(
        self, mock_read, mock_load, mock_open, invalid_config_lenghts, test_project_dir
    ):
        # Setup
        mock_load.return_value = invalid_config_lenghts
        mock_read.return_value = MagicMock(seq="ATGCATGC")

        # Execute & Verify
        with pytest.raises(ValueError) as context:
            pipeline = BasePipeline(test_project_dir)
            assert "Number of sample groups and input fastas must match!" in str(
                context.value
            )

    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_load_and_check_config_invalid_config_headers_lengths(
        self,
        mock_read,
        mock_load,
        mock_open,
        invalid_config_headers_lengths,
        test_project_dir,
    ):
        # Setup
        mock_load.return_value = invalid_config_headers_lengths
        mock_read.return_value = MagicMock(seq="ATGCATGC")

        # Execute & Verify
        with pytest.raises(ValueError) as context:
            pipeline = BasePipeline(test_project_dir)
            assert (
                "Number of header formats must be 0 or match number of sample groups!"
                in str(context.value)
            )

    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_load_and_check_config_invalid_config_headers(
        self, mock_read, mock_load, mock_open, invalid_config_headers, test_project_dir
    ):
        # Setup
        mock_load.return_value = invalid_config_headers
        mock_read.return_value = MagicMock(seq="ATGCATGC")

        # Execute & Verify
        with pytest.raises(ValueError) as context:
            pipeline = BasePipeline(test_project_dir)
            assert "Header format must be 'gisaid', 'ncbi' or 'plain'!" in str(
                context.value
            )

    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_load_and_check_config_invalid_config_mode(
        self, mock_read, mock_load, mock_open, invalid_config_mode, test_project_dir
    ):
        # Setup
        mock_load.return_value = invalid_config_mode
        mock_read.return_value = MagicMock(seq="ATGCATGC")

        # Execute & Verify
        with pytest.raises(ValueError) as context:
            pipeline = BasePipeline(test_project_dir)
            assert (
                "Unrecognized analysis mode, should by ava or ovr. Check your config file!"
                in str(context.value)
            )

    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    def test_load_and_check_config_missing_reference(
        self, mock_load, mock_open, valid_config_no_header, test_project_dir
    ):
        # Setup
        mock_load.return_value = valid_config_no_header

        # Execute & Verify
        with pytest.raises(ValueError) as context:
            pipeline = BasePipeline(test_project_dir)
            assert (
                "Missing reference sequence file or the file has incorrect format!"
                in str(context.value)
            )

    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    def test_load_and_check_config_missing_reference(
        self, mock_load, mock_open, valid_config_no_header, test_project_dir
    ):
        # Setup
        mock_load.return_value = valid_config_no_header
        with open(test_project_dir / "input/ref_seq.fasta", "w") as f:
            f.write("")

        # Execute & Verify
        with pytest.raises(ValueError) as context:
            pipeline = BasePipeline(test_project_dir)
            assert (
                "Missing reference sequence file or the file has incorrect format!"
                in str(context.value)
            )


class TestBasePipelineLoadKmerSet:
    def test_load_kmer_set_success(self, test_base_pipeline, test_project_dir):
        os.makedirs(test_project_dir / "input/indices")
        kmer_set_path = test_project_dir / "input/indices/5mer_group1_set.pkl"
        with open(kmer_set_path, "wb") as f:
            pickle.dump({"AAAAA", "CCCCC"}, f)
        kmer_set_path = test_project_dir / "input/indices/5mer_group2_set.pkl"
        with open(kmer_set_path, "wb") as f:
            pickle.dump({"GGGGG", "TTTTT"}, f)
        test_base_pipeline.load_kmer_set(test_base_pipeline.k)
        assert set(test_base_pipeline.kmer_set) == {"AAAAA", "CCCCC", "GGGGG", "TTTTT"}

    def test_load_kmer_set_missing_files(self, test_base_pipeline):
        with pytest.raises(
            FileNotFoundError,
            match="Some or all 5-mers are missing from the project directory! Please run PORTEK find_k!",
        ):
            test_base_pipeline.load_kmer_set(test_base_pipeline.k)

    def test_load_kmer_set_partial_files(self, test_base_pipeline, test_project_dir):
        os.makedirs(test_project_dir / "input/indices")
        kmer_set_path = test_project_dir / "input/indices/5mer_group1_set.pkl"
        with open(kmer_set_path, "wb") as f:
            pickle.dump({"AAAAA", "CCCCC"}, f)
        with pytest.raises(
            FileNotFoundError,
            match="Some or all 5-mers are missing from the project directory! Please run PORTEK find_k!",
        ):
            test_base_pipeline.load_kmer_set(test_base_pipeline.k)


class TestBasePipelineLoadSampleList:
    def test_load_sample_list_success(self, test_base_pipeline, test_project_dir):
        os.makedirs(test_project_dir / "input/indices")
        sample_list_path = test_project_dir / "input/indices/group1_sample_list.pkl"
        with open(sample_list_path, "wb") as f:
            pickle.dump(["sample1", "sample2"], f)
        sample_list_path = test_project_dir / "input/indices/group2_sample_list.pkl"
        with open(sample_list_path, "wb") as f:
            pickle.dump(["sample1", "sample2"], f)
        test_base_pipeline.load_sample_list()
        assert sorted(test_base_pipeline.sample_list) == [
            "group1_sample1",
            "group1_sample2",
            "group2_sample1",
            "group2_sample2",
        ]

    def test_load_sample_list_missing_files(self, test_base_pipeline):
        with pytest.raises(
            FileNotFoundError,
            match="Some or all samples are missing from the project directory! Please run PORTEK find_k!",
        ):
            test_base_pipeline.load_sample_list()

    def test_load_sample_list_partial_files(self, test_base_pipeline, test_project_dir):
        os.makedirs(test_project_dir / "input/indices")
        sample_list_path = test_project_dir / "input/indices/group1_sample_list.pkl"
        with open(sample_list_path, "wb") as f:
            pickle.dump(["sample1", "sample2"], f)
        with pytest.raises(
            FileNotFoundError,
            match="Some or all samples are missing from the project directory! Please run PORTEK find_k!",
        ):
            test_base_pipeline.load_sample_list()


class TestBasePipelineCheckMinMaxK:
    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_invalid_mink_type(
        self, mock_read, mock_load, mock_open, test_project_dir, valid_config_no_header
    ):
        # Setup
        mock_load.return_value = valid_config_no_header
        mock_read.return_value.seq = "ATGCATGC"
        pipeline = BasePipeline(test_project_dir, 5)

        # Execute & Verify
        with pytest.raises(
            TypeError, match="Minimum k must by an odd integer not smaller than 5!"
        ):
            pipeline.check_min_max_k("5", 7)

    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_mink_even(
        self, mock_read, mock_load, mock_open, test_project_dir, valid_config_no_header
    ):
        # Setup
        mock_load.return_value = valid_config_no_header
        mock_read.return_value.seq = "ATGCATGC"
        pipeline = BasePipeline(test_project_dir, 5)

        # Execute & Verify
        with pytest.raises(
            TypeError, match="Minimum k must by an odd integer not smaller than 5!"
        ):
            pipeline.check_min_max_k(mink=6, maxk=7)

    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_mink_small(
        self, mock_read, mock_load, mock_open, test_project_dir, valid_config_no_header
    ):
        # Setup
        mock_load.return_value = valid_config_no_header
        mock_read.return_value.seq = "ATGCATGC"
        pipeline = BasePipeline(test_project_dir, 5)

        # Execute & Verify
        with pytest.raises(
            TypeError, match="Minimum k must by an odd integer not smaller than 5!"
        ):
            pipeline.check_min_max_k(mink=3, maxk=7)

    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_invalid_maxk_type(
        self, mock_read, mock_load, mock_open, test_project_dir, valid_config_no_header
    ):
        # Setup
        mock_load.return_value = valid_config_no_header
        mock_read.return_value.seq = "ATGCATGC"
        pipeline = BasePipeline(test_project_dir, 5)

        # Execute & Verify
        with pytest.raises(
            TypeError, match="Maximum k must by an odd integer not smaller than 5!"
        ):
            pipeline.check_min_max_k(mink=5, maxk="7")

    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_maxk_even(
        self, mock_read, mock_load, mock_open, test_project_dir, valid_config_no_header
    ):
        # Setup
        mock_load.return_value = valid_config_no_header
        mock_read.return_value.seq = "ATGCATGC"
        pipeline = BasePipeline(test_project_dir, 5)

        # Execute & Verify
        with pytest.raises(
            TypeError, match="Maximum k must by an odd integer not smaller than 5!"
        ):
            pipeline.check_min_max_k(mink=5, maxk=8)

    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_maxk_less_than_mink(
        self, mock_read, mock_load, mock_open, test_project_dir, valid_config_no_header
    ):
        # Setup
        mock_load.return_value = valid_config_no_header
        mock_read.return_value.seq = "ATGCATGC"
        pipeline = BasePipeline(test_project_dir, 5)

        # Execute & Verify
        with pytest.raises(
            ValueError, match="Minimum k must be no greater than maximum k!"
        ):
            pipeline.check_min_max_k(mink=15, maxk=5)


class TestEncodeSeq:
    def test_encode_seq_valid_input(self):
        assert portek.encode_seq("ACGT") == ["00", "01", "10", "11"]
        assert portek.encode_seq("A") == ["00"]
        assert portek.encode_seq("C") == ["01"]
        assert portek.encode_seq("G") == ["10"]
        assert portek.encode_seq("T") == ["11"]

    def test_encode_seq_invalid_input(self):
        assert portek.encode_seq("N") == ["X"]
        assert portek.encode_seq("ACGN") == ["00", "01", "10", "X"]
        assert portek.encode_seq("") == []

    def test_encode_seq_mixed_case(self):
        assert portek.encode_seq("aCgT") == ["00", "01", "10", "11"]
