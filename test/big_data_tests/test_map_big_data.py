import pytest
import shutil
from portek.portek_map import MappingPipeline
from portek.portek_map_str import MappingPipelineString


class TestMapBigData:
    def test_map_bit_kmers(self, benchmark):
        def index_and_mapping_by_ambi():
            mapping_pipeline = MappingPipeline("projects/HIV_M_DNA", 15)
            mapping_pipeline.index_ref_seq(2, True)
            mapping_pipeline.run_mapping_by_ambi(2, True)

        benchmark(index_and_mapping_by_ambi)

    def test_map_str_kmers(self, benchmark):
        def index_and_mapping_by_ambi():
            mapping_pipeline = MappingPipelineString("projects/HIV_M_DNA", 15)
            mapping_pipeline.index_ref_seq(2, True)
            mapping_pipeline.run_mapping_by_ambi(2, True)

        benchmark(index_and_mapping_by_ambi)
