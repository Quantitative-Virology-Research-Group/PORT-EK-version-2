import pytest
import pickle
import os
from unittest.mock import patch, mock_open
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

    def test_map_kmer_to_index(
        self, mapping_pipeline: MappingPipeline, test_project_dir
    ):
        max_del_distance = 2
        mapping_pipeline.index_ref_seq(max_del_distance)
        print(mapping_pipeline.kmer_index)
        mapping_pipeline.mapping_dict = {"ACAAT": {}}
        mapping_pipeline.unmapped_counter = {d: 0 for d in range(max_del_distance + 1)}
        mapping_pipeline._map_kmer_to_index("ACAAT", 0)
        mapping_pipeline._map_kmer_to_index("ACAAT", 1)
        mapping_pipeline._map_kmer_to_index("ACAAT", 2)
        print(mapping_pipeline.mapping_dict)
        print(mapping_pipeline.unmapped_counter)
        assert False
