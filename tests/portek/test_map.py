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
    mock_read_seq.return_value.seq = "AAAAATTTTTCCCCCGGGGG"
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
            (0, "AAAAA", {"AAAAA"}),
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
            (0, "AAAAA", 0, {"AAAAA"}),
            (
                1,
                "AAAAA",
                0,
                {"ANAAAA", "AANAAA", "AAANAA", "AAAANA", "AAAAT"},
            ),
            (
                1,
                "TTCTT",
                5,
                {"TNTCTT", "TTNCTT", "TTCNTT", "TTCTNT"},
            ),
            (
                1,
                "ATTCT",
                4,
                {"ANTTCT", "ATNTCT", "ATTNCT", "ATTCNT", "ATTTT", "ATCTT"},
            ),
            (
                2,
                "ATTCT",
                4,
                {"ANNTTCT", "ATNNTCT", "ATTNNCT", "ATTCNNT"},
            ),
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

    # @pytest.mark.parametrize(
    #     "max_mismatches, expected_index",
    #     [
    #         (
    #             0,
    #             {
    #                 0: {
    #                     "AAAAA": [1],
    #                     "AAAAT": [2],
    #                     "AAATT": [3],
    #                     "AATTT": [4],
    #                     "ATTTT": [5],
    #                     "TTTTT": [6],
    #                 }
    #             },
    #         ),
    #         (
    #             1,
    #             {
    #                 0: {
    #                     "AAAAA": [1],
    #                     "AAAAT": [2],
    #                     "AAATT": [3],
    #                     "AATTT": [4],
    #                     "ATTTT": [5],
    #                     "TTTTT": [6],
    #                 },
    #                 1: {
    #                     "AAAAN": [1, 2],
    #                     "AAANA": [1],
    #                     "AANAA": [1],
    #                     "ANAAA": [1],
    #                     "NAAAA": [1],
    #                     "AAANT": [2, 3],
    #                     "AANAT": [2],
    #                     "ANAAT": [2],
    #                     "NAAAT": [2],
    #                     "AAATN": [3],
    #                     "AANTT": [3, 4],
    #                     "ANATT": [3],
    #                     "NAATT": [3],
    #                     "AATTN": [4],
    #                     "AATNT": [4],
    #                     "ANTTT": [4, 5],
    #                     "NATTT": [4],
    #                     "ATTTN": [5],
    #                     "ATTNT": [5],
    #                     "ATNTT": [5],
    #                     "ANTTT": [4, 5],
    #                     "NTTTT": [5, 6],
    #                     "TTTTN": [6],
    #                     "TTTNT": [6],
    #                     "TTNTT": [6],
    #                     "TNTTT": [6],
    #                 },
    #             },
    #         ),
    #     ],
    # )
    # # AAAAATTTTTs
    # def test_generate_index(
    #     self, mapping_pipeline: MappingPipeline, max_mismatches, expected_index
    # ):
    #     result_index = mapping_pipeline._generate_index(max_mismatches, "AAAAATTTTT")
    #     assert result_index[0] == expected_index[0]
    #     assert result_index.get(1, 0) == expected_index.get(1, 0)
