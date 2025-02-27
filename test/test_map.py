import pytest
import os
import pathlib
import pandas as pd
import unittest
from unittest.mock import patch, mock_open, MagicMock
from portek.portek_map import MappingPipeline


class TestMappingPipelineInit:
    @patch("os.path.isdir")
    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("pandas.read_csv")
    @patch("Bio.SeqIO.read")
    def test_init_success(
        self,
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

        pipeline = MappingPipeline("/fake/dir", 5)

        assert pipeline.project_dir == "/fake/dir"
        assert pipeline.k == 5
        assert pipeline.sample_groups == ["group1", "group2"]
        assert pipeline.mode == "ovr"
        assert pipeline.goi == "group1"
        assert pipeline.control_groups == ["group2"]
        assert pipeline.ref_seq_name == "reference"
        assert pipeline.ref_seq == "ATGC"
        assert pipeline.ref_genes == ["gene1", "gene2"]
        assert pipeline.avg_cols == ["group1_avg", "group2_avg"]
        assert "enriched" in pipeline.matrices
        assert pipeline.matrices["enriched"].equals(mock_enriched_csv)

    @patch("os.path.isdir")
    def test_init_invalid_directory(self, mock_isdir):
        mock_isdir.return_value = False
        with pytest.raises(NotADirectoryError):
            MappingPipeline("/fake/dir", 5)

    @patch("os.path.isdir")
    def test_init_invalid_k_type(self, mock_isdir):
        mock_isdir.return_value = True
        with pytest.raises(TypeError):
            MappingPipeline("/fake/dir", "not_an_int")

    @patch("os.path.isdir")
    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    def test_init_missing_config(self, mock_yaml_load, mock_open, mock_isdir):
        mock_isdir.return_value = True
        mock_open.side_effect = FileNotFoundError
        with pytest.raises(FileNotFoundError):
            MappingPipeline("/fake/dir", 5)

    @patch("os.path.isdir")
    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    def test_init_invalid_mode(
        self, mock_yaml_load, mock_open, mock_isdir, mock_config
    ):
        mock_isdir.return_value = True
        mock_config["mode"] = "invalid_mode"
        mock_yaml_load.return_value = mock_config
        with pytest.raises(ValueError):
            MappingPipeline("/fake/dir", 5)

    @patch("os.path.isdir")
    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("pandas.read_csv")
    @patch("Bio.SeqIO.read")
    def test_init_missing_ref_seq(
        self,
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
        mock_seqio_read.side_effect = ValueError
        mock_read_csv.return_value = mock_enriched_csv
        with pytest.raises(ValueError):
            MappingPipeline("/fake/dir", 5)

    @patch("os.path.isdir")
    @patch("builtins.open", new_callable=mock_open)
    @patch("yaml.safe_load")
    @patch("pandas.read_csv")
    @patch("Bio.SeqIO.read")
    def test_init_missing_enriched_csv(
        self,
        mock_seqio_read,
        mock_read_csv,
        mock_yaml_load,
        mock_open,
        mock_isdir,
        mock_config,
    ):
        mock_isdir.return_value = True
        mock_yaml_load.return_value = mock_config
        mock_seqio_read.return_value.seq = "ATGC"
        mock_read_csv.side_effect = FileNotFoundError
        with pytest.raises(FileNotFoundError):
            MappingPipeline("/fake/dir", 5)


class TestMappingPipelineGetSamples:
    @patch("builtins.open", new_callable=mock_open)
    @patch("pickle.load")
    @patch("pathlib.Path.glob")
    def test_get_samples_success(
        self, mock_glob, mock_pickle_load, mock_open, mock_proper_mapping_pipeline
    ):
        mock_glob.return_value = [
            pathlib.Path("/fake/dir/input/A_sample_list.pkl"),
            pathlib.Path("/fake/dir/input/B_sample_list.pkl"),
        ]
        mock_pickle_load.return_value = ["sample1", "sample2"]
        pipeline = mock_proper_mapping_pipeline
        pipeline.get_samples()

        assert pipeline.sample_list == [
            "A_sample1",
            "A_sample2",
            "B_sample1",
            "B_sample2",
        ]

    @patch("pathlib.Path.glob")
    def test_get_samples_empty(self, mock_glob, mock_proper_mapping_pipeline):
        mock_glob.return_value = []
        pipeline = mock_proper_mapping_pipeline
        with pytest.raises(FileNotFoundError):
            pipeline.get_samples()


class TestMappingPipeline_checkBowtie2Path:
    @patch("shutil.which")
    def test_check_bowtie2_path_success(self, mock_which, mock_proper_mapping_pipeline):
        mock_which.return_value = "/usr/bin/bowtie2"

        pipeline = mock_proper_mapping_pipeline
        assert pipeline._check_bowtie2_path() == "/usr/bin/bowtie2"

    @patch("shutil.which")
    def test_check_bowtie2_path_not_found(
        self, mock_which, mock_proper_mapping_pipeline
    ):
        mock_which.return_value = None
        pipeline = mock_proper_mapping_pipeline

        assert pipeline._check_bowtie2_path() == None

    
class TestMappingPipeline_checkIndexBuilt:
    @patch("pathlib.Path.glob")
    def test_check_index_built_success(self, mock_glob, mock_proper_mapping_pipeline):
        mock_glob.return_value = [pathlib.Path("/fake/dir/index.bt2")]
        pipeline = mock_proper_mapping_pipeline
        assert pipeline._check_index_built() == True

    @patch("pathlib.Path.glob")
    def test_check_index_built_not_found(self, mock_glob, mock_proper_mapping_pipeline):
        mock_glob.return_value = []
        pipeline = mock_proper_mapping_pipeline
        assert pipeline._check_index_built() == False


class TestMappingPipeline_bowtieBuildIndex:

    @patch('os.makedirs')
    @patch('os.path.exists')
    @patch('subprocess.run')
    def test_bowtie_build_index_success(self, mock_subprocess_run, mock_path_exists, mock_makedirs, mock_proper_mapping_pipeline):
        # Setup
        mock_path_exists.return_value = False
        mock_subprocess_run.return_value = MagicMock(returncode=0, stdout="Success", stderr="")
        pipeline = mock_proper_mapping_pipeline

        # Execute
        pipeline._bowtie_build_index(verbose=True)

        # Verify
        mock_path_exists.assert_called_once_with("/fake/dir/temp/ref_index/")
        mock_makedirs.assert_called_once_with("/fake/dir/temp/ref_index")
        mock_subprocess_run.assert_called_once_with(
            [
                "bowtie2-build",
                "-f",
                "/fake/dir/input/reference.fasta",
                "/fake/dir/temp/ref_index/reference",
            ],
            capture_output=True,
            text=True
        )

    @patch('os.makedirs')
    @patch('os.path.exists')
    @patch('subprocess.run')
    def test_bowtie_build_index_failed(self, mock_subprocess_run, mock_path_exists, mock_makedirs, mock_proper_mapping_pipeline):
        # Setup
        mock_path_exists.return_value = False
        mock_subprocess_run.return_value = MagicMock(returncode=1, stdout="", stderr="bowtie2 error")
        pipeline = mock_proper_mapping_pipeline

        # Execute
        with pytest.raises(Exception) as context:
            pipeline._bowtie_build_index(verbose=True)

        # Verify
        mock_path_exists.assert_called_once_with("/fake/dir/temp/ref_index/")
        mock_makedirs.assert_called_once_with("/fake/dir/temp/ref_index")
        mock_subprocess_run.assert_called_once_with(
            [
                "bowtie2-build",
                "-f",
                "/fake/dir/input/reference.fasta",
                "/fake/dir/temp/ref_index/reference",
            ],
            capture_output=True,
            text=True
        )
        assert context.exconly()=="Exception: bowtie2 error"

    @patch('os.path.exists')
    @patch('subprocess.run')
    def test_bowtie_build_index_directory_exists(self, mock_subprocess_run, mock_path_exists, mock_proper_mapping_pipeline):
        # Setup
        mock_path_exists.return_value = True
        mock_subprocess_run.return_value = MagicMock(returncode=0, stdout="Success", stderr="")
        pipeline = mock_proper_mapping_pipeline

        # Execute
        pipeline._bowtie_build_index(verbose=True)

        # Verify
        mock_path_exists.assert_called_once_with("/fake/dir/temp/ref_index/")
        mock_subprocess_run.assert_called_once_with(
            [
                "bowtie2-build",
                "-f",
                "/fake/dir/input/reference.fasta",
                "/fake/dir/temp/ref_index/reference",
            ],
            capture_output=True,
            text=True
        )


class TestMappingPipeline_bowtieMap:

    @patch('subprocess.run')
    def test_bowtie_map_success(self, mock_subprocess_run, mock_proper_mapping_pipeline):
        # Setup
        mock_subprocess_run.return_value = MagicMock(returncode=0, stdout="Success", stderr="")
        pipeline = mock_proper_mapping_pipeline

        # Execute
        pipeline._bowtie_map(verbose=True)

        # Verify
        mock_subprocess_run.assert_called_once_with(
            [
                "bowtie2",
                "--norc",
                "-a",
                "-L",
                "3",
                "-x",
                "/fake/dir/temp/ref_index/reference",
                "-f",
                "/fake/dir/temp/enriched_5mers.fasta",
                "-S",
                "/fake/dir/temp/enriched_5mers.sam",
            ],
            capture_output=True,
            text=True
        )

    @patch('subprocess.run')
    def test_bowtie_map_success(self, mock_subprocess_run, mock_proper_mapping_pipeline):
        # Setup
        mock_subprocess_run.return_value = MagicMock(returncode=1, stdout="", stderr="bowtie2 error")
        pipeline = mock_proper_mapping_pipeline

        # Execute
        with pytest.raises(Exception) as context:
            pipeline._bowtie_map(verbose=True)

        # Verify
        mock_subprocess_run.assert_called_once_with(
            [
                "bowtie2",
                "--norc",
                "-a",
                "-L",
                "3",
                "-x",
                "/fake/dir/temp/ref_index/reference",
                "-f",
                "/fake/dir/temp/enriched_5mers.fasta",
                "-S",
                "/fake/dir/temp/enriched_5mers.sam",
            ],
            capture_output=True,
            text=True
        )
        assert context.exconly() == "Exception: bowtie2 error"