import pathlib
import pickle
from unittest.mock import MagicMock, mock_open, patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

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
        mock_config_ava,
        mock_enriched_csv,
    ):
        mock_isdir.return_value = True
        mock_yaml_load.return_value = mock_config_ava
        mock_read_csv.return_value = mock_enriched_csv
        mock_seqio_read.return_value.seq = "ATGC"

        pipeline = MappingPipeline("/fake/dir", 5)

        assert pipeline.project_dir == "/fake/dir"
        assert pipeline.k == 5
        assert pipeline.sample_groups == ["group1", "group2"]
        assert pipeline.mode == "ava"
        assert pipeline.goi == None
        assert pipeline.control_groups == None
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
        self, mock_yaml_load, mock_open, mock_isdir, mock_config_ava
    ):
        mock_isdir.return_value = True
        mock_config_ava["mode"] = "invalid_mode"
        mock_yaml_load.return_value = mock_config_ava
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
        mock_config_ava,
        mock_enriched_csv,
    ):
        mock_isdir.return_value = True
        mock_yaml_load.return_value = mock_config_ava
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
        mock_config_ava,
    ):
        mock_isdir.return_value = True
        mock_yaml_load.return_value = mock_config_ava
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
                "--score-min",
                "L,-0.6,-1",
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
    def test_bowtie_map_error(self, mock_subprocess_run, mock_proper_mapping_pipeline):
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
                "--score-min",
                "L,-0.6,-1",
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


class TestMappingPipeline_runMapping:
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
                "--score-min",
                "L,-0.6,-1",
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
                "--score-min",
                "L,-0.6,-1",
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

    def test_parse_CIGAR(self, mock_proper_mapping_pipeline):
        pipeline = mock_proper_mapping_pipeline

        # Test case 1: Simple CIGAR string
        cigar_string = "10M"
        expected_output = ["M"] * 10
        assert pipeline._parse_CIGAR(cigar_string) == expected_output

        # Test case 2: CIGAR string with multiple operations
        cigar_string = "5M3I2D4M"
        expected_output = ["M"] * 5 + ["I"] * 3 + ["D"] * 2 + ["M"] * 4
        assert pipeline._parse_CIGAR(cigar_string) == expected_output

        # Test case 3: CIGAR string with all operations
        cigar_string = "2M1I1D2M"
        expected_output = ["M"] * 2 + ["I"] * 1 + ["D"] * 1 + ["M"] * 2
        assert pipeline._parse_CIGAR(cigar_string) == expected_output

        # Test case 4: Empty CIGAR string
        cigar_string = ""
        expected_output = []
        assert pipeline._parse_CIGAR(cigar_string) == expected_output

        # Test case 5: CIGAR string with large numbers
        cigar_string = "100M50I25D"
        expected_output = ["M"] * 100 + ["I"] * 50 + ["D"] * 25
        assert pipeline._parse_CIGAR(cigar_string) == expected_output


class TestMappingPipeline_readSAMToDf:

    @patch("pysam.AlignmentFile")
    def test_read_sam_to_df_success(
        self, mock_alignment_file, mock_proper_mapping_pipeline
    ):
        # Setup
        mock_read = MagicMock()
        mock_read.query_name = "AAAAA"
        mock_read.flag = 0
        mock_read.reference_start = 10
        mock_read.cigarstring = "10M"
        mock_read.get_tag.side_effect = lambda tag, _: {"NM": 1, "AS": 0}.get(tag, None)
        mock_alignment_file.return_value.__iter__.return_value = [mock_read]

        pipeline = mock_proper_mapping_pipeline

        # Execute
        df = pipeline._read_sam_to_df()

        # Verify
        expected_df = pd.DataFrame(
            {
                "kmer": ["AAAAA"],
                "flag": [0],
                "ref_pos": [11],
                "CIGAR": [["M"] * 10],
                "n_mismatch": [1],
                "score": [0],
                "group": ["group1"],
                "mutations": ["WT"],
                "mapping_ok": [1],
            }
        )
        pd.testing.assert_frame_equal(df, expected_df)

    @patch("pysam.AlignmentFile")
    def test_read_sam_to_df_ref_pos_1(
        self, mock_alignment_file, mock_proper_mapping_pipeline
    ):
        # Setup
        mock_read = MagicMock()
        mock_read.query_name = "AAAAA"
        mock_read.flag = 0
        mock_read.reference_start = 0
        mock_read.cigarstring = "10M"
        mock_read.get_tag.side_effect = lambda tag, _: {"NM": 1, "AS": 0}.get(tag, None)
        mock_alignment_file.return_value.__iter__.return_value = [mock_read]

        pipeline = mock_proper_mapping_pipeline

        # Execute
        df = pipeline._read_sam_to_df()

        # Verify
        expected_df = pd.DataFrame(
            {
                "kmer": ["AAAAA"],
                "flag": [0],
                "ref_pos": [1],
                "CIGAR": [["M"] * 10],
                "n_mismatch": [1],
                "score": [0],
                "group": ["group1"],
                "mutations": ["WT"],
                "mapping_ok": [1],
            }
        )
        pd.testing.assert_frame_equal(df, expected_df)

    @patch("pysam.AlignmentFile")
    def test_read_sam_to_df_no_aln(
        self, mock_alignment_file, mock_proper_mapping_pipeline
    ):
        # Setup
        mock_read = MagicMock()
        mock_read.query_name = "AAAAA"
        mock_read.flag = 4
        mock_read.reference_start = 0
        mock_read.cigarstring = None
        mock_read.has_tag.return_value = False
        mock_alignment_file.return_value.__iter__.return_value = [mock_read]

        pipeline = mock_proper_mapping_pipeline

        # Execute
        df = pipeline._read_sam_to_df()

        # Verify
        expected_df = pd.DataFrame(
            {
                "kmer": ["AAAAA"],
                "flag": [4],
                "ref_pos": [0],
                "CIGAR": [[]],
                "n_mismatch": [0],
                "score": [-30],
                "group": ["group1"],
                "mutations": ["WT"],
                "mapping_ok": [1],
            }
        )
        pd.testing.assert_frame_equal(df, expected_df)


class TestMappingPipeline_alingSeqs:
    def test_align_seqs_simple(self, mock_proper_mapping_pipeline):
        pipeline = mock_proper_mapping_pipeline
        ref_seq = "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAG"
        kmer = "CGTAC"
        map_pos = 3
        cigar = ["M", "M", "M", "M", "M"]

        q_seq, t_seq, ref_pos = pipeline._align_seqs(ref_seq, kmer, map_pos, cigar)

        assert q_seq == list("CGTAC")
        assert t_seq == list("CGTAC")
        assert ref_pos == [3, 4, 5, 6, 7]

    def test_align_seqs_pos1(self, mock_proper_mapping_pipeline):
        pipeline = mock_proper_mapping_pipeline
        ref_seq = "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAG"
        kmer = "ATGCG"
        map_pos = 0
        cigar = ["M", "M", "M", "M", "M"]

        q_seq, t_seq, ref_pos = pipeline._align_seqs(ref_seq, kmer, map_pos, cigar)

        assert q_seq == list("ATGCG")
        assert t_seq == list("ATGCG")
        assert ref_pos == [0, 1, 2, 3, 4]

    def test_align_seqs_with_deletion(self, mock_proper_mapping_pipeline):
        pipeline = mock_proper_mapping_pipeline
        ref_seq = "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAG"
        kmer = "CGACG"
        map_pos = 3
        cigar = ["M", "M", "D", "M", "M", "M"]

        q_seq, t_seq, ref_pos = pipeline._align_seqs(ref_seq, kmer, map_pos, cigar)

        assert q_seq == list("CG-ACG")
        assert t_seq == list("CGTACG")
        assert ref_pos == [3, 4, 5, 6, 7, 8]

    def test_align_seqs_with_insertion(self, mock_proper_mapping_pipeline):
        pipeline = mock_proper_mapping_pipeline
        ref_seq = "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAG"
        kmer = "CGTTA"
        map_pos = 3
        cigar = ["M", "M", "I", "M", "M"]

        q_seq, t_seq, ref_pos = pipeline._align_seqs(ref_seq, kmer, map_pos, cigar)

        assert q_seq == list("CGTTA")
        assert t_seq == list("CG-TA")
        assert ref_pos == [3, 4, 4, 5, 6]

    def test_align_seqs_with_complex_operations(self, mock_proper_mapping_pipeline):
        pipeline = mock_proper_mapping_pipeline
        ref_seq = "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAG"
        kmer = "CTAAC"
        map_pos = 3
        cigar = ["M", "D", "M", "I", "M", "M"]

        q_seq, t_seq, ref_pos = pipeline._align_seqs(ref_seq, kmer, map_pos, cigar)

        assert q_seq == list("C-TAAC")
        assert t_seq == list("CGT-AC")
        assert ref_pos == [3, 4, 5, 5, 6, 7]


class TestMappingPipeline_joinIndels:

    def test_join_indels_empty(self, mock_proper_mapping_pipeline):
        mutations = []
        result = mock_proper_mapping_pipeline._join_indels(mutations)
        assert result == []

    def test_join_indels_substitutions_only(self, mock_proper_mapping_pipeline):
        mutations = [(10, "A", "T"), (20, "G", "C")]
        result = mock_proper_mapping_pipeline._join_indels(mutations)
        assert result == [(10, "A", "T"), (20, "G", "C")]

    def test_join_indels_insertions_only(self, mock_proper_mapping_pipeline):
        mutations = [(10, "ins", "A"), (10, "ins", "T"), (20, "ins", "G")]
        result = mock_proper_mapping_pipeline._join_indels(mutations)
        assert result == [(10, "ins", "AT"), (20, "ins", "G")]

    def test_join_indels_deletions_only(self, mock_proper_mapping_pipeline):
        mutations = [(10, "del", 10), (11, "del", 11), (20, "del", 20)]
        result = mock_proper_mapping_pipeline._join_indels(mutations)
        assert result == [(10, "del", 11), (20, "del", 20)]

    def test_join_indels_mixed(self, mock_proper_mapping_pipeline):
        mutations = [
            (10, "A", "T"),
            (15, "ins", "A"),
            (15, "ins", "T"),
            (20, "del", 20),
            (21, "del", 21),
        ]
        result = mock_proper_mapping_pipeline._join_indels(mutations)
        assert result == [(10, "A", "T"), (15, "ins", "AT"), (20, "del", 21)]

    def test_join_indels_complex(self, mock_proper_mapping_pipeline):
        mutations = [
            (10, "A", "T"),
            (15, "ins", "A"),
            (15, "ins", "T"),
            (20, "del", 20),
            (21, "del", 21),
            (25, "ins", "G"),
            (25, "ins", "C"),
            (30, "del", 30),
            (31, "del", 31),
            (32, "G", "A"),
        ]
        result = mock_proper_mapping_pipeline._join_indels(mutations)
        assert result == [
            (10, "A", "T"),
            (15, "ins", "AT"),
            (20, "del", 21),
            (25, "ins", "GC"),
            (30, "del", 31),
            (32, "G", "A"),
        ]


class TestMappingPipeline_findVariants:
    def test_find_variants_no_mutations(self, mock_proper_mapping_pipeline):
        pipeline = mock_proper_mapping_pipeline
        ref_seq = "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAG"
        kmer = "CGTAC"
        map_pos = 4
        cigar = ["M", "M", "M", "M", "M"]

        mutations = pipeline._find_variants(ref_seq, kmer, map_pos, cigar)
        assert mutations == []

    def test_find_variants_with_substitution(self, mock_proper_mapping_pipeline):
        pipeline = mock_proper_mapping_pipeline
        ref_seq = "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAG"
        kmer = "CGTTC"
        map_pos = 4
        cigar = ["M", "M", "M", "M", "M"]

        mutations = pipeline._find_variants(ref_seq, kmer, map_pos, cigar)
        assert mutations == [(7, "A", "T")]

    def test_find_variants_with_deletion(self, mock_proper_mapping_pipeline):
        pipeline = mock_proper_mapping_pipeline
        ref_seq = "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAG"
        kmer = "CGACG"
        map_pos = 4
        cigar = ["M", "M", "D", "M", "M", "M"]

        mutations = pipeline._find_variants(ref_seq, kmer, map_pos, cigar)
        assert mutations == [(6, "del", 6)]

    def test_find_variants_with_insertion(self, mock_proper_mapping_pipeline):
        pipeline = mock_proper_mapping_pipeline
        ref_seq = "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAG"
        kmer = "CGTTA"
        map_pos = 4
        cigar = ["M", "M", "I", "M", "M"]

        mutations = pipeline._find_variants(ref_seq, kmer, map_pos, cigar)
        assert mutations == [(5, "ins", "T")]

    def test_find_variants_with_complex_operations(self, mock_proper_mapping_pipeline):
        pipeline = mock_proper_mapping_pipeline
        ref_seq = "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAG"
        kmer = "CTAAG"
        map_pos = 4
        cigar = ["M", "D", "M", "M", "M", "M"]

        mutations = pipeline._find_variants(ref_seq, kmer, map_pos, cigar)
        assert mutations == [(5, "del", 5), (8, "C", "A")]


class TestMappingPipeline_mutationTupleToText:
    def test_mutation_tuple_to_text_deletion(self, mock_proper_mapping_pipeline):
        pipeline = mock_proper_mapping_pipeline
        mutation_tuple = (10, "del", 15)
        result = pipeline._mutation_tuple_to_text(mutation_tuple)
        assert result == "10_15del"

    def test_mutation_tuple_to_text_insertion(self, mock_proper_mapping_pipeline):
        pipeline = mock_proper_mapping_pipeline
        mutation_tuple = (10, "ins", "A")
        result = pipeline._mutation_tuple_to_text(mutation_tuple)
        assert result == "10_11insA"

    def test_mutation_tuple_to_text_substitution(self, mock_proper_mapping_pipeline):
        pipeline = mock_proper_mapping_pipeline
        mutation_tuple = (10, "A", "T")
        result = pipeline._mutation_tuple_to_text(mutation_tuple)
        assert result == "10A>T"


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


class TestMappingPipeline_getNPeaks:
    def test_get_n_peaks_single_peak(self, mock_proper_mapping_pipeline):
        actual_positions = {"group1": {"AAAAA": [10]}}
        pipeline = mock_proper_mapping_pipeline
        result = pipeline._get_n_peaks("group1", "AAAAA", 10, actual_positions)
        assert result == 1

    def test_get_n_peaks_multiple_peaks(self, mock_proper_mapping_pipeline):
        actual_positions = {"group1": {"AAAAA": [10, 20, 30]}}
        pipeline = mock_proper_mapping_pipeline
        result = pipeline._get_n_peaks("group1", "AAAAA", 10, actual_positions)
        assert result == 3

    def test_get_n_peaks_no_peaks(self, mock_proper_mapping_pipeline):
        actual_positions = {"group1": {"AAAAA": []}}
        pipeline = mock_proper_mapping_pipeline
        result = pipeline._get_n_peaks("group1", "AAAAA", 10, actual_positions)
        assert result == 0


class TestMappingPipeline_resolveMultipleMappings:

    def test_resolve_multiple_mappings_single_mapping(
        self, mock_proper_mapping_pipeline
    ):
        # Setup
        data = {
            "kmer": ["AAAAA", "TTTTT"],
            "flag": [0, 0],
            "ref_pos": [10, 20],
            "score": [0, 0],
            "real_pos": [10, 20],
            "n_peaks": [1, 1],
            "mapping_ok": [1, 1],
        }
        mappings_df = pd.DataFrame(data)
        pipeline = mock_proper_mapping_pipeline

        # Execute
        bad_mappings = pipeline._resolve_multiple_mappings(mappings_df)

        # Verify
        assert bad_mappings.empty

    def test_resolve_multiple_mappings_multiple_mappings(
        self, mock_proper_mapping_pipeline
    ):
        # Setup
        data = {
            "kmer": ["AAAAA", "AAAAA", "TTTTT", "TTTTT"],
            "flag": [0, 256, 0, 256],
            "ref_pos": [10, 15, 20, 25],
            "score": [0, -6, 0, 0],
            "real_pos": [10, 10, 25, 25],
            "n_peaks": [1, 1, 1, 1],
            "mapping_ok": [1, 1, 1, 1],
        }
        mappings_df = pd.DataFrame(data)
        pipeline = mock_proper_mapping_pipeline

        # Execute
        bad_mappings = pipeline._resolve_multiple_mappings(mappings_df)

        # Verify
        assert len(bad_mappings) == 2
        assert bad_mappings.equals(pd.Index([1, 2]))

    def test_resolve_multiple_mappings_equal(self, mock_proper_mapping_pipeline):
        # Setup
        data = {
            "kmer": ["AAAAA", "AAAAA", "TTTTT", "TTTTT"],
            "flag": [0, 256, 0, 256],
            "ref_pos": [10, 15, 20, 25],
            "score": [30, 30, 40, 40],
            "real_pos": [10, 15, 20, 25],
            "n_peaks": [2, 2, 2, 2],
            "mapping_ok": [1, 1, 1, 1],
        }
        mappings_df = pd.DataFrame(data)
        pipeline = mock_proper_mapping_pipeline

        # Execute
        bad_mappings = pipeline._resolve_multiple_mappings(mappings_df)

        # Verify
        assert bad_mappings.empty


class TestMappingPipeline_filterMappings:

    @patch("portek.portek_map.MappingPipeline._resolve_multiple_mappings")
    def test_filter_mappings_bad_mappings(
        self, mock_resolve_multiple_mappings, mock_proper_mapping_pipeline
    ):
        # Setup
        pipeline = mock_proper_mapping_pipeline
        mock_resolve_multiple_mappings.return_value = pd.Index([1, 3])
        data = {
            "kmer": ["AAAAA", "AAAAA", "CCCCC", "CCCCC"],
            "flag": [0, 256, 0, 256],
            "ref_pos": [10, 20, 30, 40],
            "CIGAR": ["10M", "10M", "10M", "10M"],
            "n_mismatch": [0, 1, 0, 1],
            "score": [0, -6, 0, -6],
            "group": ["group1", "group2", "group1", "group2"],
            "mutations": ["WT", "WT", "WT", "WT"],
            "mapping_ok": [1, 1, 1, 1],
            "real_pos": [10, 10, 30, 30],
            "n_peaks": [1, 1, 1, 1],
        }
        mappings_df = pd.DataFrame(data)

        # Execute
        filtered_df = pipeline._filter_mappings(mappings_df)

        # Verify
        assert len(filtered_df) == 2
        assert filtered_df["kmer"].tolist() == ["AAAAA", "CCCCC"]

    @patch("portek.portek_map.MappingPipeline._resolve_multiple_mappings")
    def test_filter_mappings_no_bad_mappings(
        self, mock_resolve_multiple_mappings, mock_proper_mapping_pipeline
    ):
        # Setup
        pipeline = mock_proper_mapping_pipeline
        mock_resolve_multiple_mappings.return_value = pd.Index([])
        data = {
            "kmer": ["AAAAA", "AAAAA", "CCCCC", "CCCCC"],
            "flag": [0, 256, 0, 256],
            "ref_pos": [10, 20, 30, 40],
            "CIGAR": ["10M", "10M", "10M", "10M"],
            "n_mismatch": [0, 0, 0, 0],
            "score": [0, 0, 0, 0],
            "group": ["group1", "group2", "group1", "group2"],
            "mutations": ["WT", "WT", "WT", "WT"],
            "mapping_ok": [1, 1, 1, 1],
            "real_pos": [10, 20, 30, 40],
            "n_peaks": [2, 2, 2, 2],
        }
        mappings_df = pd.DataFrame(data)

        # Execute
        filtered_df = pipeline._filter_mappings(mappings_df)

        # Verify
        assert len(filtered_df) == 4
        assert filtered_df["kmer"].tolist() == ["AAAAA", "AAAAA", "CCCCC", "CCCCC"]


class TestMappingPipeline_formatMappingsDf:

    def test_format_mappings_df(self, mock_proper_mapping_pipeline):
        pipeline = mock_proper_mapping_pipeline

        # Create a sample DataFrame to test
        data = {
            "kmer": ["AAAAA", "TTTTT", "CCCCC"],
            "flag": [0, 4, 0],
            "ref_pos": [10, 0, 30],
            "CIGAR": [["M"] * 5, [], ["M"] * 5],
            "n_mismatch": [0, 0, 1],
            "score": [0, -30, -6],
            "group": ["group1", "group2", "group1"],
            "mutations": ["WT", "WT", "31G>C"],
            "mapping_ok": [1, 1, 1],
            "real_pos": [10, 20, 30],
            "n_peaks": [1, 1, 1],
        }
        mappings_df = pd.DataFrame(data)

        # Expected DataFrame after formatting
        expected_data = {
            "kmer": ["AAAAA", "TTTTT", "CCCCC"],
            "ref_pos": [10, 0, 30],
            "n_mismatch": [0, 5, 1],
            "group": ["group1", "group2", "group1"],
            "mutations": ["WT", "-", "31G>C"],
        }
        expected_df = pd.DataFrame(expected_data).sort_values("ref_pos")

        # Execute the method
        formatted_df = pipeline._format_mappings_df(mappings_df)

        # Verify
        pd.testing.assert_frame_equal(formatted_df, expected_df)


class TestMappingPipelinePredictUnmapped:

    @patch("portek.portek_map.MappingPipeline._read_sam_to_df")
    def test_predict_unmapped_success(
        self, mock_read_sam_to_df, mock_proper_mapping_pipeline
    ):
        # Setup
        mock_read_sam_to_df.return_value = pd.DataFrame(
            {
                "kmer": ["AAAAA", "AAAAC", "AAAAG", "AAAAT"] * 10,
                "flag": [0, 0, 0, 4] * 10,
                "ref_pos": [10, 20, 30, 0] * 10,
                "real_pos": [11, 21, 31, 40] * 10,
                "n_mismatch": [0, 0, 0, 0] * 10,
                "score": [0, 0, 0, -30] * 10,
                "group": ["group1", "group1", "group1", "group1"] * 10,
                "mutations": ["WT", "WT", "WT", "WT"] * 10,
                "mapping_ok": [1, 1, 1, 1] * 10,
            }
        )
        pipeline = mock_proper_mapping_pipeline

        # Execute
        result_df = pipeline._predict_unmapped(mock_read_sam_to_df.return_value)

        # Verify
        assert "pred_pos" in result_df.columns
        assert 0.0 not in result_df["pred_pos"].values
        assert "pred_err" in result_df.columns
        assert "pred_r2" in result_df.columns

    @patch("portek.portek_map.MappingPipeline._read_sam_to_df")
    def test_predict_unmapped_one_group_unmapped(
        self, mock_read_sam_to_df, mock_proper_mapping_pipeline, capsys
    ):
        # Setup
        mock_read_sam_to_df.return_value = pd.DataFrame(
            {
                "kmer": ["AAAAA", "AAAAC", "AAAAG", "AAAAT"] * 10,
                "flag": [0, 0, 0, 4] * 10,
                "ref_pos": [10, 20, 30, 0] * 10,
                "real_pos": [11, 21, 31, 40] * 10,
                "n_mismatch": [0, 0, 0, 0] * 10,
                "score": [0, 0, 0, -30] * 10,
                "group": ["group1", "group1", "group1", "group2"] * 10,
                "mutations": ["WT", "WT", "WT", "WT"] * 10,
                "mapping_ok": [1, 1, 1, 1] * 10,
            }
        )
        pipeline = mock_proper_mapping_pipeline

        # Execute
        result_df = pipeline._predict_unmapped(mock_read_sam_to_df.return_value)
        captured_output = capsys.readouterr()

        # Verify
        assert "pred_pos" in result_df.columns
        assert result_df.loc[result_df["group"] == "group2", "pred_pos"].sum() == 0
        assert "pred_err" in result_df.columns
        assert result_df.loc[result_df["group"] == "group2", "pred_err"].sum() == 0.0
        assert "pred_r2" in result_df.columns
        assert result_df.loc[result_df["group"] == "group2", "pred_r2"].sum() == 0.0
        assert (
            captured_output.out
            == "Not enough mapped k-mers for position prediction in group group2, skipping.\n"
        )

    @patch("portek.portek_map.MappingPipeline._read_sam_to_df")
    def test_predict_unmapped_no_mapped(
        self, mock_read_sam_to_df, mock_proper_mapping_pipeline
    ):
        # Setup
        mock_read_sam_to_df.return_value = pd.DataFrame(
            {
                "kmer": ["AAAAA", "AAAAC", "AAAAG", "AAAAT"] * 10,
                "flag": [4, 4, 4, 4] * 10,
                "ref_pos": [0, 0, 0, 0] * 10,
                "real_pos": [0, 0, 0, 0] * 10,
                "n_mismatch": [3, 3, 3, 3] * 10,
                "score": [-30, -30, -30, -30] * 10,
                "group": ["group1", "group1", "group1", "group1"] * 10,
                "mutations": ["WT", "WT", "WT", "WT"] * 10,
                "mapping_ok": [1, 1, 1, 1] * 10,
            }
        )
        pipeline = mock_proper_mapping_pipeline

        # Execute
        result_df = pipeline._predict_unmapped(mock_read_sam_to_df.return_value)

        # Verify
        assert "pred_pos" in result_df.columns
        assert result_df["pred_pos"].sum() == 0
        assert "pred_err" in result_df.columns
        assert result_df["pred_pos"].sum() == 0.0
        assert "pred_r2" in result_df.columns
        assert result_df["pred_r2"].sum() == 0.0


class TestMappingPipeline_countMappings:
    def test_count_mappings_all_aligned(self, mock_proper_mapping_pipeline, capsys):
        # Setup
        pipeline = mock_proper_mapping_pipeline
        mappings_df = pd.DataFrame(
            {"kmer": ["AAAAA", "AAAAC", "AAAAG", "AAAAT"], "flag": [0, 0, 0, 0]}
        )
        # Exectue
        pipeline._count_mappings(mappings_df)

        # Verify
        captured_output = capsys.readouterr()
        assert captured_output.out == (
            "\n4 out of 4 5-mers were aligned to the reference.\n"
            "Of those, 0 were aligned to more than one position.\n"
            "0 5-mers couldn't be aligned.\n"
        )

    def test_count_mappings_some_unaligned(self, mock_proper_mapping_pipeline, capsys):
        # Setup
        pipeline = mock_proper_mapping_pipeline
        mappings_df = pd.DataFrame(
            {"kmer": ["AAAAA", "AAAAC", "AAAAG", "AAAAT"], "flag": [0, 4, 0, 4]}
        )

        # Execute
        pipeline._count_mappings(mappings_df)

        # Verify
        captured_output = capsys.readouterr()
        assert captured_output.out == (
            "\n2 out of 4 5-mers were aligned to the reference.\n"
            "Of those, 0 were aligned to more than one position.\n"
            "2 5-mers couldn't be aligned.\n"
        )

    def test_count_mappings_multiple_mappings(
        self, mock_proper_mapping_pipeline, capsys
    ):
        # Setup
        pipeline = mock_proper_mapping_pipeline
        mappings_df = pd.DataFrame(
            {"kmer": ["AAAAA", "AAAAA", "AAAAG", "AAAAT"], "flag": [0, 256, 0, 0]}
        )

        # Execute
        pipeline._count_mappings(mappings_df)

        # Verify
        captured_output = capsys.readouterr()
        assert captured_output.out == (
            "\n3 out of 4 5-mers were aligned to the reference.\n"
            "Of those, 1 were aligned to more than one position.\n"
            "1 5-mers couldn't be aligned.\n"
        )

    def test_count_mappings_no_aligned(self, mock_proper_mapping_pipeline, capsys):
        # Setup
        pipeline = mock_proper_mapping_pipeline
        mappings_df = pd.DataFrame(
            {"kmer": ["AAAAA", "AAAAC", "AAAAG", "AAAAT"], "flag": [4, 4, 4, 4]}
        )

        # Execute
        pipeline._count_mappings(mappings_df)

        # Verify
        captured_output = capsys.readouterr()
        assert captured_output.out == (
            "\n0 out of 4 5-mers were aligned to the reference.\n"
            "Of those, 0 were aligned to more than one position.\n"
            "4 5-mers couldn't be aligned.\n"
        )


class TestMappingPipelineAnalyzeMapping:

    @patch("portek.portek_map.MappingPipeline._read_sam_to_df")
    @patch("portek.portek_map.MappingPipeline._get_kmer_peaks")
    @patch("portek.portek_map.MappingPipeline._filter_mappings")
    @patch("portek.portek_map.MappingPipeline._find_variants")
    @patch("portek.portek_map.MappingPipeline._mutation_tuple_to_text")
    @patch("portek.portek_map.MappingPipeline._format_mappings_df")
    @patch("portek.portek_map.MappingPipeline._count_mappings")
    def test_analyze_mapping_success(
        self,
        mock_count_mappings,
        mock_format_mappings_df,
        mock_mutation_tuple_to_text,
        mock_find_variants,
        mock_filter_mappings,
        mock_get_kmer_peaks,
        mock_read_sam_to_df,
        mock_proper_mapping_pipeline,
    ):
        # Setup
        mock_read_sam_to_df.return_value = pd.DataFrame(
            {
                "kmer": ["AAAAA", "TTTTT"],
                "flag": [0, 0],
                "ref_pos": [10, 20],
                "CIGAR": [["M"] * 5, ["M"] * 5],
                "n_mismatch": [0, 1],
                "score": [0, -6],
                "group": ["group1", "group2"],
                "mutations": ["WT", "WT"],
                "mapping_ok": [1, 1],
            }
        )
        mock_get_kmer_peaks.return_value = {
            "group1": {"AAAAA": np.array([10])},
            "group2": {"TTTTT": np.array([20])},
        }
        mock_filter_mappings.return_value = mock_read_sam_to_df.return_value
        mock_find_variants.return_value = [(21, "A", "T")]
        mock_mutation_tuple_to_text.return_value = "21A>T"
        mock_format_mappings_df.return_value = pd.DataFrame(
            {
                "kmer": ["AAAAA", "TTTTT"],
                "ref_pos": [10, 20],
                "n_mismatch": [0, 1],
                "group": ["group1", "group2"],
                "mutations": ["WT", "21A>T"],
            }
        )
        pipeline = mock_proper_mapping_pipeline

        # Execute
        pipeline.analyze_mapping(verbose=True)

        # Verify
        mock_read_sam_to_df.assert_called_once()
        mock_get_kmer_peaks.assert_called_once()
        mock_filter_mappings.assert_called_once()
        mock_find_variants.assert_called_once()
        mock_mutation_tuple_to_text.assert_called_once()
        mock_format_mappings_df.assert_called_once()
        mock_count_mappings.assert_called_once()

        assert "mappings" in pipeline.matrices
        assert pipeline.matrices["mappings"].loc[1, "mutations"] == "21A>T"

    @patch("portek.portek_map.MappingPipeline._read_sam_to_df")
    @patch("portek.portek_map.MappingPipeline._get_kmer_peaks")
    @patch("portek.portek_map.MappingPipeline._filter_mappings")
    @patch("portek.portek_map.MappingPipeline._find_variants")
    @patch("portek.portek_map.MappingPipeline._mutation_tuple_to_text")
    @patch("portek.portek_map.MappingPipeline._format_mappings_df")
    def test_analyze_mapping_no_mismatches(
        self,
        mock_format_mappings_df,
        mock_mutation_tuple_to_text,
        mock_find_variants,
        mock_filter_mappings,
        mock_get_kmer_peaks,
        mock_read_sam_to_df,
        mock_proper_mapping_pipeline,
    ):
        # Setup
        mock_read_sam_to_df.return_value = pd.DataFrame(
            {
                "kmer": ["AAAAA"],
                "flag": [0],
                "ref_pos": [10],
                "CIGAR": [["M"] * 5],
                "n_mismatch": [0],
                "score": [0],
                "group": ["group1"],
                "mutations": ["WT"],
                "mapping_ok": [1],
            }
        )
        mock_get_kmer_peaks.return_value = {"group1": {"AAAAA": np.array([10])}}
        mock_filter_mappings.return_value = mock_read_sam_to_df.return_value
        mock_format_mappings_df.return_value = pd.DataFrame(
            {
                "kmer": ["AAAAA"],
                "ref_pos": [10],
                "n_mismatch": [0],
                "group": ["group1"],
                "mutations": ["WT"],
            }
        )
        pipeline = mock_proper_mapping_pipeline

        # Execute
        pipeline.analyze_mapping(verbose=True)

        # Verify
        mock_read_sam_to_df.assert_called_once()
        mock_get_kmer_peaks.assert_called_once()
        mock_filter_mappings.assert_called_once()
        mock_find_variants.assert_not_called()
        mock_mutation_tuple_to_text.assert_not_called()
        mock_format_mappings_df.assert_called_once()

        assert "mappings" in pipeline.matrices
        assert pipeline.matrices["mappings"].loc[0, "mutations"] == "WT"

    @patch("portek.portek_map.MappingPipeline._read_sam_to_df")
    @patch("portek.portek_map.MappingPipeline._get_kmer_peaks")
    @patch("portek.portek_map.MappingPipeline._filter_mappings")
    @patch("portek.portek_map.MappingPipeline._find_variants")
    @patch("portek.portek_map.MappingPipeline._mutation_tuple_to_text")
    @patch("portek.portek_map.MappingPipeline._format_mappings_df")
    def test_analyze_mapping_with_multiple_mismatches(
        self,
        mock_format_mappings_df,
        mock_mutation_tuple_to_text,
        mock_find_variants,
        mock_filter_mappings,
        mock_get_kmer_peaks,
        mock_read_sam_to_df,
        mock_proper_mapping_pipeline,
    ):
        # Setup
        mock_read_sam_to_df.return_value = pd.DataFrame(
            {
                "kmer": ["AAAAA"],
                "flag": [0],
                "ref_pos": [10],
                "CIGAR": [["M", "M", "I", "M", "M"]],
                "n_mismatch": [2],
                "score": [-12],
                "group": ["group1"],
                "mutations": ["WT"],
                "mapping_ok": [1],
            }
        )
        mock_get_kmer_peaks.return_value = {"group1": {"AAAAA": np.array([10])}}
        mock_filter_mappings.return_value = mock_read_sam_to_df.return_value
        mock_find_variants.return_value = [(10, "T", "A"), (11, "ins", "A")]
        mock_mutation_tuple_to_text.side_effect = ["10T>A", "11_12insA"]
        mock_format_mappings_df.return_value = pd.DataFrame(
            {
                "kmer": ["AAAAA"],
                "ref_pos": [10],
                "n_mismatch": [2],
                "group": ["group1"],
                "mutations": ["10T>A; 11_12insA"],
            }
        )
        pipeline = mock_proper_mapping_pipeline

        # Execute
        pipeline.analyze_mapping(verbose=True)

        # Verify
        mock_read_sam_to_df.assert_called_once()
        mock_get_kmer_peaks.assert_called_once()
        mock_filter_mappings.assert_called_once()
        mock_find_variants.assert_called_once()
        mock_mutation_tuple_to_text.assert_called()
        mock_format_mappings_df.assert_called_once()

        assert "mappings" in pipeline.matrices
        assert pipeline.matrices["mappings"].loc[0, "mutations"] == "10T>A; 11_12insA"


class TestMappingPipelineSaveMappingsDf:

    @patch("pandas.DataFrame.to_csv")
    def test_save_mappings_df_success(self, mock_to_csv, mock_proper_mapping_pipeline):
        # Setup
        pipeline = mock_proper_mapping_pipeline
        pipeline.matrices["mappings"] = pd.DataFrame(
            {
                "kmer": ["AAAAA", "TTTTT"],
                "ref_pos": [1, 2],
                "n_mismatch": [0, 0],
                "group": ["group1", "group2"],
                "mutations": ["WT", "WT"],
                "pred_pos": [2, 3],
                "pred_err": [1, 1],
                "pred_r2": [0.99, 0.99],
            }
        )

        # Execute
        pipeline.save_mappings_df()

        # Verify
        mock_to_csv.assert_called_once_with(
            f"{pipeline.project_dir}/output/enriched_{pipeline.k}mers_mappings.csv"
        )

    @patch("pandas.DataFrame.to_csv")
    def test_save_mappings_df_no_mappings(
        self, mock_to_csv, mock_proper_mapping_pipeline
    ):
        # Setup
        pipeline = mock_proper_mapping_pipeline
        pipeline.matrices["mappings"] = pd.DataFrame()

        # Execute
        pipeline.save_mappings_df()

        # Verify
        mock_to_csv.assert_called_once_with(
            f"{pipeline.project_dir}/output/enriched_{pipeline.k}mers_mappings.csv"
        )


class TestMappingPipeline_setHistogramAxProperties:

    @patch("matplotlib.axes.Axes.set_title")
    @patch("matplotlib.axes.Axes.set_xlim")
    @patch("matplotlib.axes.Axes.set_xlabel")
    @patch("matplotlib.axes.Axes.set_ylabel")
    def test_set_ax_properites_first_in_row(
        self,
        mock_set_ylabel,
        mock_set_xlabel,
        mock_set_xlim,
        mock_set_title,
        mock_proper_mapping_pipeline,
    ):
        # Setup
        pipeline = mock_proper_mapping_pipeline
        _, ax = plt.subplots()

        # Exectue
        pipeline._set_histogram_ax_properties(ax, 1000, "group1", 0, 3)

        # Verify
        assert mock_set_title.call_args[0][0] == "group1"
        assert mock_set_xlim.call_args[0] == (0, 1000)
        assert mock_set_xlabel.call_args[0][0] == "Reference genome position"
        assert mock_set_ylabel.called
        assert mock_set_ylabel.call_args[0][0] == "K-mer counts"

    @patch("matplotlib.axes.Axes.set_title")
    @patch("matplotlib.axes.Axes.set_xlim")
    @patch("matplotlib.axes.Axes.set_xlabel")
    @patch("matplotlib.axes.Axes.set_ylabel")
    def test_set_ax_properites_not_first_in_row(
        self,
        mock_set_ylabel,
        mock_set_xlabel,
        mock_set_xlim,
        mock_set_title,
        mock_proper_mapping_pipeline,
    ):
        # Setup
        pipeline = mock_proper_mapping_pipeline
        _, ax = plt.subplots()

        # Exectue
        pipeline._set_histogram_ax_properties(ax, 1000, "group1", 1, 3)

        # Verify
        assert mock_set_title.call_args[0][0] == "group1"
        assert mock_set_xlim.call_args[0] == (0, 1000)
        assert mock_set_xlabel.call_args[0][0] == "Reference genome position"
        assert not mock_set_ylabel.called

    @patch("matplotlib.axes.Axes.set_title")
    @patch("matplotlib.axes.Axes.set_xlim")
    @patch("matplotlib.axes.Axes.set_xlabel")
    @patch("matplotlib.axes.Axes.set_ylabel")
    def test_set_ax_properites_only(
        self,
        mock_set_ylabel,
        mock_set_xlabel,
        mock_set_xlim,
        mock_set_title,
        mock_proper_mapping_pipeline,
    ):
        # Setup
        pipeline = mock_proper_mapping_pipeline
        _, ax = plt.subplots()

        # Exectue
        pipeline._set_histogram_ax_properties(ax, 1000, "group1", 0, 1)

        # Verify
        assert mock_set_title.call_args[0][0] == "group1"
        assert mock_set_xlim.call_args[0] == (0, 1000)
        assert mock_set_xlabel.call_args[0][0] == "Reference genome position"
        assert mock_set_ylabel.called
        assert mock_set_ylabel.call_args[0][0] == "K-mer counts"


class TestMappingPipelinePlotKmerHistograms:

    @patch("portek.portek_map.MappingPipeline._set_histogram_ax_properties")
    @patch("matplotlib.pyplot.savefig")
    def test_plot_kmer_histograms(
        self, mock_savefig, mock_set_histogram_ax_properites, mock_proper_mapping_pipeline
    ):
        pipeline = mock_proper_mapping_pipeline
        pipeline.matrices["mappings"] = pd.DataFrame(
            {
                "kmer": ["AAAAA", "TTTTT"],
                "ref_pos": [1, 2],
                "n_mismatch": [0, 0],
                "group": ["group1", "group2"],
                "mutations": ["WT", "WT"],
                "pred_pos": [2, 3],
                "pred_err": [1, 1],
                "pred_r2": [0.99, 0.99],
            }
        )
        # Execute
        pipeline.plot_kmer_histograms()

        # Verify
        assert mock_savefig.called
        assert mock_savefig.call_args[0][0] == "/fake/dir/output/enriched_5-mers_coverage_histograms.svg"
        assert mock_savefig.call_args[1]["format"] == "svg"
        assert mock_savefig.call_args[1]["dpi"] == 300
        assert mock_savefig.call_args[1]["bbox_inches"] == "tight"

    @patch("seaborn.histplot")
    @patch("matplotlib.pyplot.subplots")
    @patch("portek.portek_map.MappingPipeline._set_histogram_ax_properties")
    @patch("matplotlib.pyplot.savefig")
    def test_plot_kmer_histograms_1row_1col(
        self,
        mock_savefig,
        mock_set_histogram_ax_properites,
        mock_subplots,
        mock_histplot,
        mock_proper_mapping_pipeline,
    ):
        pipeline = mock_proper_mapping_pipeline
        pipeline.matrices["mappings"] = pd.DataFrame(
            {
                "kmer": ["AAAAA", "TTTTT"],
                "ref_pos": [1, 2],
                "n_mismatch": [0, 0],
                "group": ["group1", "group1"],
                "mutations": ["WT", "WT"],
                "pred_pos": [2, 3],
                "pred_err": [1, 1],
                "pred_r2": [0.99, 0.99],
            }
        )
        # Mock subplots to return a figure and axes
        fig = MagicMock()
        ax = MagicMock()
        mock_subplots.return_value = (fig, np.array([[ax]]))

        # Execute
        pipeline.plot_kmer_histograms()

        # Verify
        assert mock_histplot.called
        assert mock_subplots.called
        assert mock_subplots.call_args[0][0] == 1
        assert mock_subplots.call_args[0][1] == 1
        assert mock_subplots.call_args[1]["figsize"] == (6, 6)
        assert mock_savefig.called

    @patch("seaborn.histplot")
    @patch("matplotlib.pyplot.subplots")
    @patch("portek.portek_map.MappingPipeline._set_histogram_ax_properties")
    @patch("matplotlib.pyplot.savefig")
    def test_plot_kmer_histograms_1row_maxcol(
        self,
        mock_savefig,
        mock_set_histogram_ax_properites,
        mock_subplots,
        mock_histplot,
        mock_proper_mapping_pipeline,
    ):
        pipeline = mock_proper_mapping_pipeline
        pipeline.matrices["mappings"] = pd.DataFrame(
            {
                "kmer": ["AAAAA", "TTTTT", "CCCCC"],
                "ref_pos": [1, 2, 3],
                "n_mismatch": [0, 0, 0],
                "group": ["group1", "group2", "group3"],
                "mutations": ["WT", "WT", "WT"],
                "pred_pos": [2, 3, 4],
                "pred_err": [1, 1, 1],
                "pred_r2": [0.99, 0.99, 0.99],
            }
        )
        # Mock subplots to return a figure and axes
        fig = MagicMock()
        ax = MagicMock()
        mock_subplots.return_value = (fig, np.array([[ax, ax, ax]]))

        # Execute
        pipeline.plot_kmer_histograms()

        # Verify
        assert mock_histplot.called
        assert mock_subplots.called
        assert mock_subplots.call_args[0][0] == 1
        assert mock_subplots.call_args[0][1] == 3
        assert mock_subplots.call_args[1]["figsize"] == (18, 6)
        assert mock_savefig.called

    @patch("seaborn.histplot")
    @patch("matplotlib.pyplot.subplots")
    @patch("portek.portek_map.MappingPipeline._set_histogram_ax_properties")
    @patch("matplotlib.pyplot.savefig")
    def test_plot_kmer_histograms_2rows_maxcol(
        self,
        mock_savefig,
        mock_set_histogram_ax_properites,
        mock_subplots,
        mock_histplot,
        mock_proper_mapping_pipeline,
    ):
        pipeline = mock_proper_mapping_pipeline
        pipeline.matrices["mappings"] = pd.DataFrame(
            {
                "kmer": ["AAAAA", "TTTTT", "CCCCC", "GGGGG"],
                "ref_pos": [1, 2, 3, 4],
                "n_mismatch": [0, 0, 0, 0],
                "group": ["group1", "group2", "group3", "group4"],
                "mutations": ["WT", "WT", "WT", "WT"],
                "pred_pos": [2, 3, 4, 5],
                "pred_err": [1, 1, 1, 1],
                "pred_r2": [0.99, 0.99, 0.99, 0.99],
            }
        )
        # Mock subplots to return a figure and axes
        fig = MagicMock()
        ax = MagicMock()
        mock_subplots.return_value = (fig, np.array([[ax, ax, ax], [ax, ax, ax]]))

        # Execute
        pipeline.plot_kmer_histograms()

        # Verify
        assert mock_histplot.called
        assert mock_subplots.called
        assert mock_subplots.call_args[0][0] == 2
        assert mock_subplots.call_args[0][1] == 3
        assert mock_subplots.call_args[1]["figsize"] == (18, 12)
        assert mock_savefig.called
