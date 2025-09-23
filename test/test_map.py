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
        "n,bit_kmer",
        [
            (0, 69905),
            (1, 69906),
            (2, 69922),
            (3, 70178),
            (4, 74274),
            (5, 139810),
            (6, 139816),
            (7, 139912),
            (8, 141448),
            (9, 166024),
            (10, 559240),
            (11, 559236),
            (12, 559172),
            (13, 558148),
            (14, 541764),
            (15, 279620),
        ],
    )
    def test_get_bit_kmer(self, mapping_pipeline: MappingPipeline, n, bit_kmer):
        bit_ref_seq = portek.encode_seq_as_bits(mapping_pipeline.ref_seq)
        result_kmer = mapping_pipeline._get_bit_kmer(bit_ref_seq, n)
        assert result_kmer == bit_kmer

    @pytest.mark.parametrize(
        "del_distance,bit_kmer,expected_ambi_kmers",
        [
            (0, 0b00010001000100010001, {0b00010001000100010001}),
            (
                1,
                0b00010001000100010001,
                {
                    0b11110001000100010001,
                    0b00011111000100010001,
                    0b00010001111100010001,
                    0b00010001000111110001,
                    0b00010001000100011111,
                },
            ),
            (
                2,
                0b00010001000100010001,
                {
                    0b11111111000100010001,
                    0b11110001111100010001,
                    0b11110001000111110001,
                    0b11110001000100011111,
                    0b00011111111100010001,
                    0b00011111000111110001,
                    0b00011111000100011111,
                    0b00010001111111110001,
                    0b00010001111100011111,
                    0b00010001000111111111,
                },
            ),
        ],
    )
    def test_generate_ambi_kmers(
        self,
        mapping_pipeline: MappingPipeline,
        del_distance,
        bit_kmer,
        expected_ambi_kmers,
    ):
        result = mapping_pipeline._generate_ambi_kmers(del_distance, bit_kmer)
        assert result == expected_ambi_kmers

    @pytest.mark.parametrize(
        "max_mismatches, expected_index",
        [
            (
                0,
                {
                    0: {
                        69905: [1],
                        69906: [2],
                        69922: [3],
                        70178: [4],
                        74274: [5],
                        139810: [6],
                    }
                },
            ),
            (
                1,
                {
                    0: {
                        69905: [1],
                        69906: [2],
                        69922: [3],
                        70178: [4],
                        74274: [5],
                        139810: [6],
                    },
                    1: {
                        0b00010001000100011111: [1, 2],
                        0b00010001000111110001: [1],
                        0b00010001111100010001: [1],
                        0b00011111000100010001: [1],
                        0b11110001000100010001: [1],
                        0b00010001000111110010: [2, 3],
                        0b00010001111100010010: [2],
                        0b00011111000100010010: [2],
                        0b11110001000100010010: [2],
                        0b00010001000100101111: [3],
                        0b00010001111100100010: [3, 4],
                        0b00011111000100100010: [3],
                        0b11110001000100100010: [3],
                        0b00010001001000101111: [4],
                        0b00010001001011110010: [4],
                        0b00011111001000100010: [4, 5],
                        0b11110001001000100010: [4],
                        0b00010010001000101111: [5],
                        0b00010010001011110010: [5],
                        0b00010010111100100010: [5],
                        0b00011111001000100010: [4, 5],
                        0b11110010001000100010: [5, 6],
                        0b00100010001000101111: [6],
                        0b00100010001011110010: [6],
                        0b00100010111100100010: [6],
                        0b00101111001000100010: [6],
                    },
                },
            ),
        ],
    )
    # AAAAATTTTTs
    def test_generate_index(
        self, mapping_pipeline: MappingPipeline, max_mismatches, expected_index
    ):
        bit_ref_seq = portek.encode_seq_as_bits(
            "AAAAATTTTT"
        )  # shorter ref_seq for smaller index
        result_index = mapping_pipeline._generate_index(max_mismatches, bit_ref_seq)
        assert result_index[0] == expected_index[0]
        assert result_index.get(1, 0) == expected_index.get(1, 0)

    def test_map_kmer_to_index_by_ambi(
        self, mapping_pipeline: MappingPipeline, benchmark
    ):
        mapping_pipeline.index_ref_seq(2)
        mapping_pipeline.mapping_dict = {
            kmer: {} for kmer in mapping_pipeline.matrices["enriched"].index
        }

        def run_mapping():
            mapping_pipeline._map_kmer_to_index_by_ambi("AAAAA", 0)
            mapping_pipeline._map_kmer_to_index_by_ambi("AAAAT", 0)
            mapping_pipeline._map_kmer_to_index_by_ambi("TTTTT", 0)
            mapping_pipeline._map_kmer_to_index_by_ambi("AAAAA", 1)
            mapping_pipeline._map_kmer_to_index_by_ambi("TTTCC", 2)

        benchmark(run_mapping)

        assert mapping_pipeline.mapping_dict["AAAAA"][0] == {1}
        assert mapping_pipeline.mapping_dict["AAAAT"][0] == {2}
        assert mapping_pipeline.mapping_dict["TTTTT"][0] == {6}
        assert mapping_pipeline.mapping_dict["AAAAA"][1] == {1, 2}
        assert mapping_pipeline.mapping_dict["TTTCC"][2] == {6, 7, 8, 9, 10}

    def test_map_kmer_to_index_bitwise(
        self, mapping_pipeline: MappingPipeline, benchmark
    ):
        mapping_pipeline.index_ref_seq(2)
        mapping_pipeline.mapping_dict = {
            kmer: {} for kmer in mapping_pipeline.matrices["enriched"].index
        }

        def run_mapping():
            mapping_pipeline._map_kmer_to_index_bitwise("AAAAA", 0)
            mapping_pipeline._map_kmer_to_index_bitwise("AAAAT", 0)
            mapping_pipeline._map_kmer_to_index_bitwise("TTTTT", 0)
            mapping_pipeline._map_kmer_to_index_bitwise("AAAAA", 1)
            mapping_pipeline._map_kmer_to_index_bitwise("TTTCC", 2)

        benchmark(run_mapping)

        assert mapping_pipeline.mapping_dict["AAAAA"][0] == {1}
        assert mapping_pipeline.mapping_dict["AAAAT"][0] == {2}
        assert mapping_pipeline.mapping_dict["TTTTT"][0] == {6}
        assert mapping_pipeline.mapping_dict["AAAAA"][1] == {1, 2}
        assert mapping_pipeline.mapping_dict["TTTCC"][2] == {6, 7, 8, 9, 10}
