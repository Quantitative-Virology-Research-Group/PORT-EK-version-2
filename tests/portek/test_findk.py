import pytest
import os
import yaml
from unittest.mock import patch, mock_open
from portek.portek_findk import KmerFinder, FindOptimalKPipeline

class TestKmerFinderInit:
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_valid_init(self, mock_read, mock_load, test_project_dir, valid_config_no_header):
        #Setup
        mock_load.return_value = valid_config_no_header
        mock_read.return_value.seq = "ATGCATGC"
        os.makedirs(test_project_dir / "input")
        with open(test_project_dir / "config.yaml", "w") as f:
            yaml.safe_dump(valid_config_no_header, f)
        with open(test_project_dir / "input/group1.fasta", "w") as f:
            f.write(">seq1\nATGCATGC\n>seq2\nATGCATGC")
        with open(test_project_dir / "input/group2.fasta", "w") as f:
            f.write(">seq1\nATGCATGC\n>seq2\nATGCATGC")

        #Execute
        kmer_finder = KmerFinder(test_project_dir, mink=5, maxk=7)

        #Verify
        assert kmer_finder.mink == 5
        assert kmer_finder.maxk == 7
        assert len(kmer_finder.seq_lists) == 2
        for seq_list in kmer_finder.seq_lists:
            assert [seq.id for seq in seq_list] == ["seq1", "seq2"]

    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_valid_init_headers(self, mock_read, mock_load, test_project_dir, valid_config):
        #Setup
        mock_load.return_value = valid_config
        mock_read.return_value.seq = "ATGCATGC"
        os.makedirs(test_project_dir / "input")
        with open(test_project_dir / "config.yaml", "w") as f:
            yaml.safe_dump(valid_config, f)
        with open(test_project_dir / "input/group1.fasta", "w") as f:
            f.write(">des/cription|seq1|date\nATGCATGC\n>des/cription|seq2|date\nATGCATGC")
        with open(test_project_dir / "input/group2.fasta", "w") as f:
            f.write(">seq1_|description\nATGCATGC\n>seq2_|description\nATGCATGC")

        #Execute
        kmer_finder = KmerFinder(test_project_dir, mink=5, maxk=7)

        #Verify
        assert kmer_finder.mink == 5
        assert kmer_finder.maxk == 7
        assert len(kmer_finder.seq_lists) == 2
        for seq_list in kmer_finder.seq_lists:
            assert [seq.id for seq in seq_list] == ["seq1", "seq2"]

    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_invalid_init_headers(self, mock_read, mock_load, test_project_dir, valid_config_no_header):
        #Setup
        mock_load.return_value = valid_config_no_header
        mock_read.return_value.seq = "ATGCATGC"
        os.makedirs(test_project_dir / "input")
        with open(test_project_dir / "config.yaml", "w") as f:
            yaml.safe_dump(valid_config_no_header, f)
        with open(test_project_dir / "input/group1.fasta", "w") as f:
            f.write(">des/cription|seq1|date\nATGCATGC\n>des/cription|seq2|date\nATGCATGC")
        with open(test_project_dir / "input/group2.fasta", "w") as f:
            f.write(">seq1_|description\nATGCATGC\n>seq2_|description\nATGCATGC")

        #Execute & Verify
        with pytest.raises(ValueError, match="Sequence ids cannot contain '/'. Specify correct header formats in the config file."):
            KmerFinder(test_project_dir, mink=5, maxk=7)


class TestFindOptimalKPipeline:
    @patch("yaml.safe_load")
    @patch("Bio.SeqIO.read")
    def test_valid_init(self, mock_read, mock_load, test_project_dir, valid_config_no_header):
        #Setup
        mock_load.return_value = valid_config_no_header
        mock_read.return_value.seq = "ATGCATGC"
        times = {5:10}
        with open(test_project_dir / "config.yaml", "w") as f:
            yaml.safe_dump(valid_config_no_header, f)

        #Execute
        findk_pipeline = FindOptimalKPipeline(test_project_dir, mink=5, maxk=7, times=times)

        #Verify
        assert findk_pipeline.mink == 5
        assert findk_pipeline.maxk == 7
        assert findk_pipeline.times == {5:10}
        assert findk_pipeline.avg_cols == ["group1_avg","group2_avg"]
        assert findk_pipeline.freq_cols == ["group1_freq","group2_freq"]
        assert findk_pipeline.c_cols == ["group1_c","group2_c"]
        assert findk_pipeline.f_cols == ["group1_f","group2_f"]