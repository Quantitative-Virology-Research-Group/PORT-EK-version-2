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

    def test_map_kmer_to_index(
        self, mapping_pipeline: MappingPipeline, test_project_dir
    ):
        max_del_distance = 2
        mapping_pipeline.index_ref_seq(max_del_distance)
        mapping_pipeline.mapping_dict = {"ACAAT": {}}
        mapping_pipeline.unmapped_counter = {d: 0 for d in range(max_del_distance + 1)}
        mapping_pipeline._map_kmer_to_index("ACAAT", 0)
        mapping_pipeline._map_kmer_to_index("ACAAT", 1)
        mapping_pipeline._map_kmer_to_index("ACAAT", 2)
