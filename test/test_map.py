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
        seldf,
        mapping_pipeline: MappingPipeline,
        del_distance,
        bit_kmer,
        expected_ambi_kmers,
    ):
        result = mapping_pipeline._generate_ambi_kmers(del_distance, bit_kmer)
        assert result == expected_ambi_kmers
