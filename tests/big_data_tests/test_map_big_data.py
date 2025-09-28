import pytest
import shutil
from portek.portek_map import MappingPipeline


class TestMapBigData:
    def test_map_bit_kmers(self, benchmark):
        def index_and_mapping_by_ambi():
            mapping_pipeline = MappingPipeline("projects/HIV_M_DNA", 15)
            mapping_pipeline.index_ref_seq(2, True)
            mapping_pipeline.run_mapping(2, True)

        benchmark(index_and_mapping_by_ambi)
