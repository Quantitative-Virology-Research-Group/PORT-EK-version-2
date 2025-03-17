import pytest
import pandas as pd
import numpy as np
from portek import EnrichedKmersPipeline


class TestEnrichedKmersPipeline:

    def test_compare_group_pair(self, mock_enriched_pipeline_ava):
        # Setup
        pipeline = mock_enriched_pipeline_ava
        pipeline.sample_group_dict = {
            "group1": ["group1_sample1", "group1_sample2"],
            "group2": ["group2_sample1", "group2_sample2"],
        }
        pipeline.matrices = {
            "test_matrix": pd.DataFrame(
                {
                    "group1_sample1": [1, 1, 0],
                    "group1_sample2": [0, 1, 0],
                    "group2_sample1": [1, 0, 1],
                    "group2_sample2": [0, 0, 1],
                    "group1_avg": [0.5, 1, 0],
                    "group2_avg": [0.5, 0, 1],
                },
                index=["kmer1", "kmer2", "kmer3"],
            )
        }
        matrix_type = "test_matrix"
        group1 = "group1"
        group2 = "group2"
        verbose = False

        result = pipeline._compare_group_pair(matrix_type, group1, group2, verbose)

        assert result[0] == group1
        assert result[1] == group2
        assert isinstance(result[2], pd.Series)
        assert isinstance(result[3], pd.Series)
        assert isinstance(result[4], pd.Series)
        assert result[2].equals(
            pipeline.matrices[matrix_type]["group1_avg"]
            - pipeline.matrices[matrix_type]["group2_avg"]
        )
        assert len(result[3]) == len(pipeline.matrices[matrix_type])
        assert len(result[4]) == len(pipeline.matrices[matrix_type])
