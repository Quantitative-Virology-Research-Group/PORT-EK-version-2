import pytest
import os
import pathlib
import pickle
import pandas as pd
import numpy as np
import unittest
from unittest.mock import patch, mock_open, MagicMock
from scipy.stats import linregress
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

    @patch("os.makedirs")
    @patch("os.path.exists")
    @patch("subprocess.run")
    def test_bowtie_build_index_success(
        self,
        mock_subprocess_run,
        mock_path_exists,
        mock_makedirs,
        mock_proper_mapping_pipeline,
    ):
        # Setup
        mock_path_exists.return_value = False
        mock_subprocess_run.return_value = MagicMock(
            returncode=0, stdout="Success", stderr=""
        )
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
            text=True,
        )

    @patch("os.makedirs")
    @patch("os.path.exists")
    @patch("subprocess.run")
    def test_bowtie_build_index_failed(
        self,
        mock_subprocess_run,
        mock_path_exists,
        mock_makedirs,
        mock_proper_mapping_pipeline,
    ):
        # Setup
        mock_path_exists.return_value = False
        mock_subprocess_run.return_value = MagicMock(
            returncode=1, stdout="", stderr="bowtie2 error"
        )
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
            text=True,
        )
        assert context.exconly() == "Exception: bowtie2 error"

    @patch("os.path.exists")
    @patch("subprocess.run")
    def test_bowtie_build_index_directory_exists(
        self, mock_subprocess_run, mock_path_exists, mock_proper_mapping_pipeline
    ):
        # Setup
        mock_path_exists.return_value = True
        mock_subprocess_run.return_value = MagicMock(
            returncode=0, stdout="Success", stderr=""
        )
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
            text=True,
        )


class TestMappingPipeline_bowtieMap:

    @patch("subprocess.run")
    def test_bowtie_map_success(
        self, mock_subprocess_run, mock_proper_mapping_pipeline
    ):
        # Setup
        mock_subprocess_run.return_value = MagicMock(
            returncode=0, stdout="Success", stderr=""
        )
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
            text=True,
        )

    @patch("subprocess.run")
    def test_bowtie_map_success(
        self, mock_subprocess_run, mock_proper_mapping_pipeline
    ):
        # Setup
        mock_subprocess_run.return_value = MagicMock(
            returncode=1, stdout="", stderr="bowtie2 error"
        )
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
            text=True,
        )
        assert context.exconly() == "Exception: bowtie2 error"


class TestMappingPipelineRunMapping:
    @patch("portek.portek_map.MappingPipeline._check_bowtie2_path")
    @patch("portek.portek_map.MappingPipeline._check_index_built")
    @patch("subprocess.run")
    def test_run_mapping_success(
        self,
        mock_subprocess_run,
        mock_check_index_built,
        mock_check_bowtie2_path,
        mock_proper_mapping_pipeline,
    ):
        # Setup
        mock_check_bowtie2_path.return_value = "/usr/bin/bowtie2"
        mock_check_index_built.return_value = True
        mock_subprocess_run.return_value = MagicMock(
            returncode=0, stdout="Success", stderr=""
        )
        pipeline = mock_proper_mapping_pipeline

        # Execute
        pipeline.run_mapping()

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
            text=True,
        )

    @patch("portek.portek_map.MappingPipeline._check_bowtie2_path")
    def test_run_mapping_no_bowtie(
        self,
        mock_check_bowtie2_path,
        mock_proper_mapping_pipeline,
    ):
        # Setup
        mock_check_bowtie2_path.return_value = None
        pipeline = mock_proper_mapping_pipeline

        # Execute
        with pytest.raises(FileNotFoundError) as context:
            pipeline.run_mapping()

        # Verify
        assert (
            context.exconly()
            == "FileNotFoundError: bowtie2 not found! Please install bowtie and add it to your PATH!"
        )

    @patch("portek.portek_map.MappingPipeline._check_bowtie2_path")
    @patch("portek.portek_map.MappingPipeline._check_index_built")
    @patch("os.makedirs")
    @patch("subprocess.run")
    def test_run_mapping_index_not_exists(
        self,
        mock_subprocess_run,
        mock_makedirs,
        mock_check_index_built,
        mock_check_bowtie2_path,
        mock_proper_mapping_pipeline,
    ):
        # Setup
        mock_check_bowtie2_path.return_value = "/usr/bin/bowtie2"
        mock_check_index_built.return_value = False

        mock_subprocess_run.return_value = MagicMock(
            returncode=0, stdout="Success", stderr=""
        )
        pipeline = mock_proper_mapping_pipeline

        # Exectue
        pipeline.run_mapping()

        # Verify
        mock_makedirs.assert_called_once_with("/fake/dir/temp/ref_index")
        mock_subprocess_run.assert_any_call(
            [
                "bowtie2-build",
                "-f",
                "/fake/dir/input/reference.fasta",
                "/fake/dir/temp/ref_index/reference",
            ],
            capture_output=True,
            text=True,
        )
        mock_subprocess_run.assert_called_with(
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
            text=True,
        )


class TestMappingPipeline_parseCIGAR:

    def test_parse_CIGAR_success(
        self,
    ):
        assert False

    def test_parse_CIGAR_empty(
        self,
    ):
        assert False


class TestMappingPipeline_readSAMToDf:

    def test_read_sam_to_df_success(
        self,
    ):
        assert False

    def test_read_sam_to_df_no_file(
        self,
    ):
        assert False

    def test_read_sam_to_df_wrong_format(
        self,
    ):
        assert False


class TestMappingPipeline_loadKmerPos:

    @patch("builtins.open", new_callable=mock_open)
    @patch("pickle.load")
    def test_load_kmer_pos_success(
        self, mock_pickle_load, mock_open, mock_proper_mapping_pipeline
    ):
        # Setup
        mock_pickle_load.return_value = {0: [1, 2, 3], 1: [4, 5, 6]}
        pipeline = mock_proper_mapping_pipeline
        kmer_distros = {}
        # Execute
        pipeline._load_kmer_pos(
            pathlib.Path("/fake/dir/input/indices/5mers_group1_pos_dict.pkl"),
            ["AAAAA", "AAAAC"],
            kmer_distros,
        )

        # Verify
        mock_open.assert_called_once_with(
            pathlib.Path("/fake/dir/input/indices/5mers_group1_pos_dict.pkl"), mode="rb"
        )
        assert kmer_distros == {
            "group1_enriched": {"AAAAA": [1, 2, 3], "AAAAC": [4, 5, 6]}
        }

    @patch("builtins.open", new_callable=mock_open)
    def test_load_kmer_pos_file_not_found(
        self, mock_open, mock_proper_mapping_pipeline
    ):
        # Setup
        mock_open.side_effect = FileNotFoundError
        pipeline = mock_proper_mapping_pipeline
        kmer_distros = {}

        # Execute & Verify
        with pytest.raises(FileNotFoundError):
            pipeline._load_kmer_pos(
                pathlib.Path("/fake/dir/input/indices/5mers_group1_pos_dict.pkl"),
                ["AAAAA", "AAAAC"],
                kmer_distros,
            )

    @patch("builtins.open", new_callable=mock_open)
    @patch("pickle.load")
    def test_load_kmer_pos_empty_file(
        self, mock_pickle_load, mock_open, mock_proper_mapping_pipeline
    ):
        # Setup
        mock_pickle_load.return_value = {}
        pipeline = mock_proper_mapping_pipeline
        kmer_distros = {}

        # Execute
        pipeline._load_kmer_pos(
            pathlib.Path("/fake/dir/input/indices/5mers_group1_pos_dict.pkl"),
            ["AAAAA", "AAAAC"],
            kmer_distros,
        )

        # Verify
        mock_open.assert_called_once_with(
            pathlib.Path("/fake/dir/input/indices/5mers_group1_pos_dict.pkl"), mode="rb"
        )
        assert kmer_distros == {"group1_enriched": {}}

    @patch("builtins.open", new_callable=mock_open)
    @patch("pickle.load")
    def test_load_kmer_pos_invalid_format(
        self, mock_pickle_load, mock_open, mock_proper_mapping_pipeline
    ):
        # Setup
        mock_pickle_load.side_effect = pickle.UnpicklingError
        pipeline = mock_proper_mapping_pipeline
        kmer_distros = {}

        # Execute & Verify
        with pytest.raises(pickle.UnpicklingError):
            pipeline._load_kmer_pos(
                pathlib.Path("/fake/dir/input/indices/5mers_group1_pos_dict.pkl"),
                ["AAAAA", "AAAAC"],
                kmer_distros,
            )


class TestMappingPipeline_getTotalDistros:

    def test_get_total_distros_success(self, mock_proper_mapping_pipeline):
        # Setup
        pipeline = mock_proper_mapping_pipeline
        kmer_distros = {"group1_enriched": {"AAAAA": [1, 2, 3], "AAAAC": [4, 5, 6]}}

        # Execute
        pipeline._get_total_distros(kmer_distros)

        # Verify
        assert kmer_distros == {
            "group1_enriched": {"AAAAA": [1, 2, 3], "AAAAC": [4, 5, 6]},
            "conserved": {"AAAAA": [1, 2, 3], "AAAAC": [4, 5, 6]},
        }

    def test_get_total_distros_empty(self, mock_proper_mapping_pipeline):
        # Setup
        pipeline = mock_proper_mapping_pipeline
        kmer_distros = {}

        # Execute
        pipeline._get_total_distros(kmer_distros)

        # Verify
        assert kmer_distros == {"conserved": {}}

    def test_get_total_distros_multiple_groups(self, mock_proper_mapping_pipeline):
        # Setup
        pipeline = mock_proper_mapping_pipeline
        kmer_distros = {
            "group1_enriched": {"AAAAA": [1, 2, 3], "AAAAC": [4, 5, 6]},
            "group2_enriched": {"AAAAG": [7, 8, 9], "AAAAT": [10, 11, 12]},
        }

        # Execute
        pipeline._get_total_distros(kmer_distros)

        # Verify
        assert kmer_distros == {
            "group1_enriched": {"AAAAA": [1, 2, 3], "AAAAC": [4, 5, 6]},
            "group2_enriched": {"AAAAG": [7, 8, 9], "AAAAT": [10, 11, 12]},
            "conserved": {
                "AAAAA": [1, 2, 3],
                "AAAAC": [4, 5, 6],
                "AAAAG": [7, 8, 9],
                "AAAAT": [10, 11, 12],
            },
        }

    def test_get_total_distros_overlapping_kmers(self, mock_proper_mapping_pipeline):
        # Setup
        pipeline = mock_proper_mapping_pipeline
        kmer_distros = {
            "group1_enriched": {"AAAAA": [1, 2, 3], "AAAAC": [4, 5, 6]},
            "group2_enriched": {"AAAAA": [7, 8, 9], "AAAAT": [10, 11, 12]},
        }

        # Execute
        pipeline._get_total_distros(kmer_distros)

        # Verify
        assert kmer_distros == {
            "group1_enriched": {"AAAAA": [1, 2, 3], "AAAAC": [4, 5, 6]},
            "group2_enriched": {"AAAAA": [7, 8, 9], "AAAAT": [10, 11, 12]},
            "conserved": {
                "AAAAA": [1, 2, 3, 7, 8, 9],
                "AAAAC": [4, 5, 6],
                "AAAAT": [10, 11, 12],
            },
        }


class TestMappingPipeline_getPeaksFromDistro:

    def test_get_peaks_from_distro_one_peak(self, mock_proper_mapping_pipeline):
        # Setup
        pipeline = mock_proper_mapping_pipeline
        kmer_distros = {
            "group1": {
                "AAAAA": [
                    1,
                    2,
                    2,
                    3,
                    4,
                    5,
                    6,
                    7,
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    15,
                    16,
                    17,
                    18,
                    19,
                    20,
                ]
            }
        }
        expected_peaks = np.array([2])

        # Execute
        peaks = pipeline._get_peaks_from_distro(kmer_distros)

        # Verify
        assert (peaks["group1"]["AAAAA"] == expected_peaks).all()

    def test_get_peaks_from_distro_two_close_peaks(self, mock_proper_mapping_pipeline):
        # Setup
        pipeline = mock_proper_mapping_pipeline
        kmer_distros = {
            "group1": {
                "AAAAA": [
                    1,
                    2,
                    2,
                    2,
                    3,
                    4,
                    4,
                    5,
                    6,
                    7,
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    15,
                    16,
                    17,
                    18,
                    19,
                    20,
                ]
            }
        }
        expected_peaks = np.array([2])

        # Execute
        peaks = pipeline._get_peaks_from_distro(kmer_distros)

        # Verify
        assert (peaks["group1"]["AAAAA"] == expected_peaks).all()

    def test_get_peaks_from_distro_two_peaks(self, mock_proper_mapping_pipeline):
        # Setup
        pipeline = mock_proper_mapping_pipeline
        kmer_distros = {
            "group1": {
                "AAAAA": [
                    1,
                    2,
                    2,
                    2,
                    3,
                    4,
                    5,
                    6,
                    7,
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    15,
                    16,
                    16,
                    16,
                    17,
                    18,
                    19,
                    20,
                ]
            }
        }
        expected_peaks = np.array([2, 16])

        # Execute
        peaks = pipeline._get_peaks_from_distro(kmer_distros)

        # Verify
        assert (peaks["group1"]["AAAAA"] == expected_peaks).all()

    def test_get_peaks_from_distro_empty(self, mock_proper_mapping_pipeline):
        # Setup
        pipeline = mock_proper_mapping_pipeline
        kmer_distros = {"conserved": {}}
        expected_peaks = {"conserved": {}}

        # Execute
        peaks = pipeline._get_peaks_from_distro(kmer_distros)

        # Verify
        assert peaks == expected_peaks

    def test_get_peaks_from_distro_no_peaks(self, mock_proper_mapping_pipeline):
        # Setup
        pipeline = mock_proper_mapping_pipeline
        kmer_distros = {"group1": {"AAAAA": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]}}
        expected_peaks = np.array([])

        # Execute
        peaks = pipeline._get_peaks_from_distro(kmer_distros)

        # Verify
        assert (peaks["group1"]["AAAAA"] == expected_peaks).all()


class TestMappingPipeline_getKmerPeaks:
    def test_get_kmer_peaks(self, mock_proper_mapping_pipeline):
        # Setup
        data = {
            "kmer": ["AAAAA", "TTTTT", "CCCCC"],
            "group": ["group1", "group2", "group1"],
        }
        mappings_df = pd.DataFrame(data)
        mapping_pipeline = mock_proper_mapping_pipeline

        # Mock dependencies
        mock_files = [MagicMock(spec=pathlib.Path) for _ in range(2)]
        for i, mock_file in enumerate(mock_files):
            mock_file.stem = f"{i}_pos_dict"
        with patch.object(pathlib.Path, "glob", return_value=mock_files), patch.object(
            MappingPipeline, "_load_kmer_pos"
        ) as mock_load_kmer_pos, patch.object(
            MappingPipeline, "_get_total_distros"
        ) as mock_get_total_distros, patch.object(
            MappingPipeline,
            "_get_peaks_from_distro",
            return_value={
                "group1": {"AAAAA": [1, 2], "CCCCC": [3, 4]},
                "group2": {"TTTTT": [5, 6]},
            },
        ) as mock_get_peaks_from_distro:
            # Execute
            distro_peaks = mapping_pipeline._get_kmer_peaks(mappings_df, verbose=True)

            # Assertions
            mock_load_kmer_pos.assert_called()
            mock_get_total_distros.assert_called_once()
            mock_get_peaks_from_distro.assert_called_once()
            assert distro_peaks == {
                "group1": {"AAAAA": [1, 2], "CCCCC": [3, 4]},
                "group2": {"TTTTT": [5, 6]},
            }


class TestMappingPipeline_getRealPos:
    def test_get_real_pos_exact_match(
        self, mock_proper_mapping_pipeline, actual_positions
    ):
        pipeline = mock_proper_mapping_pipeline
        result = pipeline._get_real_pos("group1", "AAAAA", 20, actual_positions)
        assert result == 20

    def test_get_real_pos_closest_match(
        self, mock_proper_mapping_pipeline, actual_positions
    ):
        pipeline = mock_proper_mapping_pipeline
        result = pipeline._get_real_pos("group1", "AAAAA", 22, actual_positions)
        assert result == 20

    def test_get_real_pos_first_peak(
        self, mock_proper_mapping_pipeline, actual_positions
    ):
        pipeline = mock_proper_mapping_pipeline
        result = pipeline._get_real_pos("group2", "GGGGG", 6, actual_positions)
        assert result == 5

    def test_get_real_pos_last_peak(
        self, mock_proper_mapping_pipeline, actual_positions
    ):
        pipeline = mock_proper_mapping_pipeline
        result = pipeline._get_real_pos("group2", "TTTTT", 29, actual_positions)
        assert result == 30

    def test_get_real_pos_middle_peak(
        self, mock_proper_mapping_pipeline, actual_positions
    ):
        pipeline = mock_proper_mapping_pipeline
        result = pipeline._get_real_pos("group1", "CCCCC", 26, actual_positions)
        assert result == 25

    def test_get_real_pos_multiple_groups(
        self, mock_proper_mapping_pipeline, actual_positions
    ):
        pipeline = mock_proper_mapping_pipeline
        result1 = pipeline._get_real_pos("group1", "AAAAA", 20, actual_positions)
        result2 = pipeline._get_real_pos("group2", "TTTTT", 20, actual_positions)
        assert result1 == 20
        assert result2 == 20


class TestMappinPipeline_predictPos:
    def test_predict_pos(self, mock_proper_mapping_pipeline):
        pipeline = mock_proper_mapping_pipeline

        ref_pos = pd.Series([11, 12, 13, 14, 15])
        real_pos = pd.Series([1, 2, 3, 4, 5])
        regress = linregress(real_pos, ref_pos)

        pred_pos, pred_err = pipeline._predict_pos(real_pos, ref_pos, regress)

        expected_pred_pos = pd.Series([11, 12, 13, 14, 15], name="pred_pos").astype(
            float
        )
        expected_pred_err = pd.Series([0, 0, 0, 0, 0], name="pred_err").astype(float)

        pd.testing.assert_series_equal(pred_pos, expected_pred_pos)
        pd.testing.assert_series_equal(pred_err, expected_pred_err)


class TestMappingPipeline_tuneRegress:
    def test_tune_regress_success(self, mock_proper_mapping_pipeline):
        ref_pos = pd.Series(
            [11, 12, 13, 14, 15, 511, 512, 513, 514, 515, 1011, 1012, 53, 504, 1515]
        )
        real_pos = pd.Series(
            [1, 2, 3, 4, 5, 501, 502, 503, 504, 505, 1001, 1002, 1003, 1004, 1005]
        )
        pipeline = mock_proper_mapping_pipeline
        regress, thr = pipeline._tune_regress(real_pos, ref_pos)
        assert thr == 300
        assert regress.slope == 1.0
        assert regress.intercept == 10.0

    def test_tune_regress_thr100(self, mock_proper_mapping_pipeline):
        ref_pos = pd.Series(
            [11, 12, 13, 14, 15, 511, 512, 513, 514, 515, 1011, 1012, 1213, 814, 1215]
        )
        real_pos = pd.Series(
            [1, 2, 3, 4, 5, 501, 502, 503, 504, 505, 1001, 1002, 1003, 1004, 1005]
        )
        pipeline = mock_proper_mapping_pipeline
        regress, thr = pipeline._tune_regress(real_pos, ref_pos)
        assert thr == 100
        assert regress.slope == 1.0
        assert regress.intercept == 10.0

    def test_tune_regress_thrmax(self, mock_proper_mapping_pipeline):
        ref_pos = pd.Series(
            [11, 12, 13, 14, 15, 511, 512, 513, 514, 515, 1011, 1012, 1013, 1014, 2015]
        )
        real_pos = pd.Series(
            [1, 2, 3, 4, 5, 501, 502, 503, 504, 505, 1001, 1002, 1003, 1004, 1005]
        )
        pipeline = mock_proper_mapping_pipeline
        regress, thr = pipeline._tune_regress(real_pos, ref_pos)
        assert thr == 800
        assert regress.slope == 1.0
        assert regress.intercept == 10.0

    def test_tune_regress_exactmatch(self, mock_proper_mapping_pipeline):
        ref_pos = pd.Series([11, 12, 13, 14, 15, 1011, 1012, 1013, 1014, 1015])
        real_pos = pd.Series([1, 2, 3, 4, 5, 1001, 1002, 1003, 1004, 1005])
        pipeline = mock_proper_mapping_pipeline
        regress, thr = pipeline._tune_regress(real_pos, ref_pos)
        assert thr == 100
        assert regress.slope == 1.0
        assert regress.intercept == 10.0
