import pytest
import pickle
import os
from unittest.mock import patch, mock_open
import portek
from portek.portek_map import MappingPipeline


@pytest.fixture()
@patch("builtins.open", new_callable=mock_open)
@patch("pandas.read_csv")
@patch("yaml.safe_load")
@patch("Bio.SeqIO.read")
def mapping_pipeline(
    mock_read_seq,
    mock_yaml_load,
    mock_read_csv,
    mock_open,
    test_project_dir,
    valid_config,
    mock_enriched_kmer_stats_csv,
) -> MappingPipeline:
    mock_yaml_load.return_value = valid_config
    mock_read_csv.return_value = mock_enriched_kmer_stats_csv
    mock_read_seq.return_value.seq = "ATCGAAA"
    return MappingPipeline(test_project_dir, 5)


class TestMappingPipelineIndexing:
    def test_check_index_exists(
        self, mapping_pipeline: MappingPipeline, test_project_dir
    ):
        os.makedirs(test_project_dir / "temp" / "ref_index")
        with open(
            f"{test_project_dir}/temp/ref_index/{"ref_seq"}_index_{5}_{2}.pkl", "wb"
        ) as f:
            pickle.dump(True, f)
        assert mapping_pipeline._check_index(2)

    def test_check_index_does_not_exist(self, mapping_pipeline: MappingPipeline):
        assert not mapping_pipeline._check_index(2)

    @pytest.mark.parametrize(
        "n_mismatches,kmer,expected_ambi_kmers",
        [
            (
                1,
                "AAAAA",
                {
                    "NAAAA",
                    "ANAAA",
                    "AANAA",
                    "AAANA",
                    "AAAAN",
                },
            ),
            (
                2,
                "AAAAA",
                {
                    "NNAAA",
                    "NANAA",
                    "NAANA",
                    "NAAAN",
                    "ANNAA",
                    "ANANA",
                    "ANAAN",
                    "AANNA",
                    "AANAN",
                    "AAANN",
                },
            ),
        ],
    )
    def test_generate_sub_kmers(
        self,
        mapping_pipeline: MappingPipeline,
        n_mismatches,
        kmer,
        expected_ambi_kmers,
    ):
        result = mapping_pipeline._generate_sub_kmers(n_mismatches, kmer)
        assert result == expected_ambi_kmers

    @pytest.mark.parametrize(
        "n_mismatches,kmer,kmer_start_position,expected_ambi_kmers",
        [
            (
                1,
                "AAAAA",
                0,
                {"ANAAA", "AANAA", "AAANA", "AAAAT"},
            ),
            (
                1,
                "TTCTT",
                5,
                {"TNTCT", "TTNCT", "TTCNT"},
            ),
            (
                1,
                "ATTCT",
                4,
                {"ANTTC", "ATNTC", "ATTNC", "ATTTT", "ATCTT"},
            ),
            (
                2,
                "ATTCT",
                4,
                {"ANNTT", "ATNNT"},
            ),
            (2, "TTCTT", 5, {"TNNTC", "TTNNC"}),
        ],
    )
    def test_generate_indel_kmers(
        self,
        mapping_pipeline: MappingPipeline,
        n_mismatches,
        kmer,
        kmer_start_position,
        expected_ambi_kmers,
    ):
        result = mapping_pipeline._generate_indel_kmers(
            n_mismatches, kmer, kmer_start_position, "AAAAATTCTT"
        )
        assert result == expected_ambi_kmers

    @pytest.mark.parametrize(
        "max_mismatches,expected_kmer_index",
        [
            (
                0,
                {
                    0: {
                        "ATCGA": [1],
                        "TCGAA": [2],
                        "CGAAA": [3],
                    }
                },
            ),
            (
                1,
                {
                    0: {
                        "ATCGA": [1],
                        "TCGAA": [2],
                        "CGAAA": [3],
                    },
                    1: {
                        "NTCGA": [1],
                        "ANCGA": [1],
                        "ATNGA": [1],
                        "ATCNA": [1],
                        "ATCGN": [1],
                        "ACGAA": [1],
                        "ATGAA": [1],
                        "ATCAA": [1],
                        "ANTCG": [1],
                        "ATNCG": [1],
                        "ATCNG": [1],
                        "NCGAA": [2],
                        "TNGAA": [2],
                        "TCNAA": [2],
                        "TCGNA": [2],
                        "TCGAN": [2],
                        "TGAAA": [2],
                        "TCAAA": [2],
                        "TCGAA": [2],
                        "TNCGA": [2],
                        "TCNGA": [2],
                        "TCGNA": [2],
                        "NGAAA": [3],
                        "CNAAA": [3],
                        "CGNAA": [3],
                        "CGANA": [3],
                        "CGAAN": [3],
                        "CNGAA": [3],
                        "CGNAA": [3],
                        "CGANA": [3],
                    },
                },
            ),
        ],
    )
    # ATCGAAA
    def test_generate_index(
        self, mapping_pipeline: MappingPipeline, max_mismatches, expected_kmer_index
    ):
        kmer_index = mapping_pipeline._generate_index(
            max_mismatches, mapping_pipeline.ref_seq
        )
        assert kmer_index == expected_kmer_index


class TestMappingPipelineMapping:
    @pytest.mark.parametrize(
        "kmer,n_mismatches,expected_dict",
        [
            ("AAAAA", 0, {0: set(), 1: set()}),
            ("ATCGA", 0, {0: {1}, 1: set()}),
            ("ATCGA", 1, {0: set(), 1: {1}}),
            ("TCGAA", 1, {0: set(), 1: {2}}),
            ("TCGAA", 0, {0: {2}, 1: set()}),
        ],
    )
    def test_map_kmer_to_index(
        self, mapping_pipeline: MappingPipeline, kmer, n_mismatches, expected_dict
    ):
        mapping_pipeline._generate_and_write_index(1, verbose=False)
        mapping_dict = {n: set() for n in range(2)}
        mapping_pipeline._map_kmer_to_index(kmer, n_mismatches, mapping_dict)
        assert mapping_dict == expected_dict

    def test_run_mapping(self, mapping_pipeline: MappingPipeline):
        mapping_pipeline._generate_and_write_index(1, verbose=False)
        print(mapping_pipeline.kmer_index)
        mapping_pipeline.run_mapping(1, verbose=False)
        assert "reference_sequence_position" in mapping_pipeline.mapping_df.columns
        assert (
            mapping_pipeline.mapping_df.at["TTCGA", "reference_sequence_position"]
            == "1,2"
        )
        assert (
            mapping_pipeline.mapping_df.at["TACGA", "reference_sequence_position"]
            == "2"
        )
        assert (
            mapping_pipeline.mapping_df.at["CGAAA", "reference_sequence_position"]
            == "3"
        )
