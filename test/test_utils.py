import pytest
import pandas as pd
from scipy import stats
from portek import calc_kmer_pvalue

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