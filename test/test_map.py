import pathlib
import pickle
import pandas as pd
from unittest.mock import patch, mock_open, MagicMock
from portek.portek_map import MappingPipeline
import pytest


class TestPortekMap:
    def test_get_kmer_pos(self):
        sample_groups = ["group1", "group2"]
        enriched_kmers = ["AAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAC"]
        sample_list = ["group1_sample1", "group1_sample2", "group2_sample1", "group2_sample2"]
        sample_group_dict = {
            "group1": ["group1_sample1", "group1_sample2"],
            "group2": ["group2_sample1", "group2_sample2"]
        }
        matrices = {
            "enriched": pd.DataFrame(index=enriched_kmers)
        }
        mapping_pipeline = MappingPipeline("test/testproject/", k=15)
        mapping_pipeline.sample_groups = sample_groups
        mapping_pipeline.sample_list = sample_list
        mapping_pipeline.sample_group_dict = sample_group_dict
        mapping_pipeline.matrices = matrices

        mock_kmer_pos_dict = {
            0: [1, 2],
            1: [4, 5]
        }


        with patch("pathlib.Path.glob", return_value=[pathlib.Path("/mock/path/15mer_group1_pos_dict.pkl"), pathlib.Path("/mock/path/15mer_group2_pos_dict.pkl")]), \
            patch("builtins.open", mock_open(read_data=pickle.dumps(mock_kmer_pos_dict))):
            result = mapping_pipeline._get_kmer_pos(verbose=False)

        expected_columns = ["group1_avg_pos", "group2_avg_pos", "total_avg_pos"]
        assert all(col in result.columns for col in expected_columns)
        assert result.shape[0] == len(enriched_kmers)
        assert result.loc["AAAAAAAAAAAAAAA", "group1_avg_pos"] == 1.5
        assert result.loc["AAAAAAAAAAAAAAA", "group2_avg_pos"] == 1.5
        assert result.loc["AAAAAAAAAAAAAAA", "total_avg_pos"] == 1.5
        assert result.loc["AAAAAAAAAAAAAAC", "group1_avg_pos"] == 4.5
        assert result.loc["AAAAAAAAAAAAAAC", "group2_avg_pos"] == 4.5
        assert result.loc["AAAAAAAAAAAAAAC", "total_avg_pos"] == 4.5

    def test_match_mappings_happy_case(self):
        ref_pos = pd.Series([0, 1000, 1500, 2000, 3000, 4000, 5000, 6000])
        avg_pos = pd.Series([100, 1100, 1800, 2100, 3100, 4100, 5100, 6700])
        expected_thr = 100
        expected_slope = 1.0
        expected_intercept = -100.0
        expected_rvalue = 1.0

        mapping_pipeline = MappingPipeline("test/testproject/", k=15)
        thr, slope, intercept, rvalue = mapping_pipeline._match_mappings(
            ref_pos, avg_pos
        )

        assert thr == expected_thr
        assert pytest.approx(slope, 0.01) == expected_slope
        assert pytest.approx(intercept, 0.01) == expected_intercept
        assert pytest.approx(rvalue, 0.01) == expected_rvalue

    def test_match_mappings_perfect_match(self):
        ref_pos = pd.Series([0, 1000, 2000, 3000, 4000, 5000])
        avg_pos = pd.Series([100, 1100, 2100, 3100, 4100, 5100])
        expected_thr = 100
        expected_slope = 1.0
        expected_intercept = -100.0
        expected_rvalue = 1.0

        mapping_pipeline = MappingPipeline("test/testproject/", k=15)
        thr, slope, intercept, rvalue = mapping_pipeline._match_mappings(
            ref_pos, avg_pos
        )

        assert thr == expected_thr
        assert pytest.approx(slope, 0.01) == expected_slope
        assert pytest.approx(intercept, 0.01) == expected_intercept
        assert pytest.approx(rvalue, 0.01) == expected_rvalue

    def test_match_mappings_no_match(self):
        ref_pos = pd.Series([0, 0, 0, 0, 4000, 5000], name="ref_pos")
        avg_pos = pd.Series([100, 1100, 2100, 3100, 4100, 5100])
        expected_thr = 0
        expected_slope = 0.0
        expected_intercept = 0.0
        expected_rvalue = 0.0

        mapping_pipeline = MappingPipeline("test/testproject/", k=15)
        thr, slope, intercept, rvalue = mapping_pipeline._match_mappings(
            ref_pos, avg_pos
        )

        assert thr == expected_thr
        assert pytest.approx(slope, 0.01) == expected_slope
        assert pytest.approx(intercept, 0.01) == expected_intercept
        assert pytest.approx(rvalue, 0.01) == expected_rvalue

    @patch.object(MappingPipeline, "_read_sam_to_df")
    @patch.object(MappingPipeline, "_get_kmer_pos")
    def test_analyze_mapping(self, mock_get_kmer_pos, mock_read_sam_to_df):
        # Mock the return values of the methods
        mock_read_sam_to_df.return_value = pd.DataFrame(
            {
                "kmer": ["AAA", "CCC", "GGG"],
                "flag": [0, 0, 0],
                "ref_pos": [100, 200, 300],
                "CIGAR": ["3M", "3M", "3M"],
                "n_mismatch": [0, 1, 2],
                "group": ["group1_enriched", "group2_enriched", "group1_enriched"],
                "mutations": ["WT", "WT", "WT"],
                "ref_pos_pred": [0.0, 0.0, 0.0],
                "ref_pos_err": [
                    0.0,
                    0.0,
                    0.0,
                ],
                "mapping_ok": [0,0,0]
            }
        ).set_index("kmer")

        mock_get_kmer_pos.return_value = pd.DataFrame(
            {
                "group1_avg_pos": [150, 250, 350],
                "group2_avg_pos": [160, 260, 360],
                "total_avg_pos": [155, 255, 355],
            },
            index=["AAA", "CCC", "GGG"],
        )

        # Call the method
        mapping_pipeline = MappingPipeline("test/testproject/", k=15)
        mapping_pipeline.analyze_mapping(verbose=True)

        # Check if the matrices["mappings"] is set correctly
        assert "mappings" in mapping_pipeline.matrices
        mappings_df = mapping_pipeline.matrices["mappings"]
        assert not mappings_df.empty
        assert "mutations" in mappings_df.columns
        assert mappings_df["mutations"].tolist() == ["WT", "WT", "WT"]
        assert mappings_df["ref_pos_pred"].tolist() == [0.0, 0.0, 0.0]
