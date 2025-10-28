import pytest
import pickle
import os
from pathlib import Path
from unittest.mock import patch, mock_open

import pandas as pd

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
    test_enriched_kmer_stats_csv,
) -> MappingPipeline:
    mock_yaml_load.return_value = valid_config
    mock_read_csv.return_value = test_enriched_kmer_stats_csv
    mock_read_seq.return_value.seq = "ATCGAAA"
    return MappingPipeline(test_project_dir, 5)


class TestMappingPipelineInit:
    def test_init(
        self,
        mapping_pipeline: MappingPipeline,
        valid_config,
        test_project_dir,
        test_enriched_kmer_stats_csv,
    ):
        assert mapping_pipeline.project_dir == test_project_dir
        assert mapping_pipeline.ref_seq == "ATCGAAA"
        assert mapping_pipeline.k == 5
        assert mapping_pipeline.sample_groups == valid_config["sample_groups"]
        assert mapping_pipeline.goi == valid_config["goi"]
        assert mapping_pipeline.ref_genes == valid_config["ref_genes"]
        assert mapping_pipeline.matrices["enriched"].index.equals(
            test_enriched_kmer_stats_csv.index
        )
        assert mapping_pipeline.kmer_index == {}
        assert mapping_pipeline.matrices["mapping"].columns.tolist() == [
            "reference_sequence_position",
            "gene",
            "number_of_mismatches",
            "group",
            "exclusivity",
        ]
        assert mapping_pipeline.matrices["mapping"].index.equals(
            mapping_pipeline.matrices["enriched"].index
        )
        assert mapping_pipeline.matrices["coverage"].columns.tolist() == [
            "group1_enriched_kmer_coverage",
            "group2_enriched_kmer_coverage",
            "conserved_kmer_coverage",
        ]
        assert mapping_pipeline.matrices["coverage"].index.equals(
            pd.Index([1, 2, 3, 4, 5, 6, 7])
        )
        assert mapping_pipeline.matrices["coverage"].sum().sum() == 0


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
    def test_generate_index(
        self, mapping_pipeline: MappingPipeline, max_mismatches, expected_kmer_index
    ):
        kmer_index = mapping_pipeline._generate_index(
            max_mismatches, mapping_pipeline.ref_seq  # type: ignore
        )
        assert kmer_index == expected_kmer_index

    def test_write_index(self, test_project_dir, mapping_pipeline: MappingPipeline):
        test_kmer_index = {0: {"AAAAA": [1]}, 1: {"AANAA": [1]}}
        max_mismatches = 1

        mapping_pipeline._write_index(max_mismatches, test_kmer_index)

        assert os.path.exists(f"{test_project_dir}/temp/ref_index")

        index_path = f"{test_project_dir}/temp/ref_index/ref_seq_index_5_1.pkl"
        assert os.path.exists(index_path)
        with open(index_path, "rb") as f:
            loaded_index = pickle.load(f)
        assert loaded_index == test_kmer_index

    def test_generate_and_write_index(
        self, mapping_pipeline: MappingPipeline, test_project_dir
    ):
        mapping_pipeline._generate_and_write_index(1, verbose=False)
        index_path = f"{test_project_dir}/temp/ref_index/ref_seq_index_5_1.pkl"
        assert os.path.exists(index_path)
        with open(index_path, "rb") as f:
            loaded_index = pickle.load(f)
        expected_index = {
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
        }
        assert loaded_index == expected_index

    def test_generate_and_write_index_verbosity(
        self, mapping_pipeline: MappingPipeline, capsys
    ):
        mapping_pipeline._generate_and_write_index(1, verbose=True)
        captured = capsys.readouterr()
        assert (
            "No 5-mer index with maximum distance 1 exists for ref_seq, building index."
            in captured.out
        )
        assert (
            "Finished building 5-mer index with maximum distance 1 exists for ref_seq"
            in captured.out
        )

    def test_load_index(self, mapping_pipeline: MappingPipeline):
        test_kmer_index = {0: {"AAAAA": [1]}, 1: {"AANAA": [1]}}
        max_mismatches = 1
        mapping_pipeline._write_index(max_mismatches, test_kmer_index)
        mapping_pipeline._load_index(max_mismatches, verbose=False)
        assert mapping_pipeline.kmer_index == test_kmer_index

    def test_load_index_verbosity(self, mapping_pipeline: MappingPipeline, capsys):
        test_kmer_index = {0: {"AAAAA": [1]}, 1: {"AANAA": [1]}}
        max_mismatches = 1
        mapping_pipeline._write_index(max_mismatches, test_kmer_index)
        mapping_pipeline._load_index(max_mismatches, verbose=True)
        captured = capsys.readouterr()
        assert (
            "5-mer index with maximum distance 1 for ref_seq already exists, loading..."
            in captured.out
        )


class TestMappingPipelineMapping:
    @pytest.fixture(scope="class")
    def expected_mapping_df(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "reference_sequence_position": ["1,2", "2", "3"],
                "gene": ["gene1, gene2", "gene2", ""],
                "number_of_mismatches": [1, 1, 0],
                "group": [
                    "group2_enriched",
                    "group1_enriched",
                    "conserved",
                ],
                "exclusivity": [
                    "exclusive",
                    "exclusive",
                    "non-exclusive",
                ],
            },
            index=["TTCGA", "TACGA", "CGAAA"],
            columns=[
                "reference_sequence_position",
                "gene",
                "number_of_mismatches",
                "group",
                "exclusivity",
            ],
        )

    @pytest.fixture(scope="class")
    def expected_coverage_df(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "group1_enriched_kmer_coverage": [0, 1, 1, 1, 1, 1, 0],
                "group2_enriched_kmer_coverage": [1, 2, 2, 2, 2, 1, 0],
                "conserved_kmer_coverage": [0, 0, 1, 1, 1, 1, 1],
            },
            index=pd.Index(
                [pos for pos in range(1, 8)], name="reference_sequence_position"
            ),
        )

    @pytest.mark.parametrize(
        "kmer,n_mismatches,expected_dict",
        [
            ("AAAAA", 0, {0: set()}),
            ("ATCGA", 0, {0: {1}}),
            ("ATCGA", 1, {0: {1}, 1: {1}}),
            ("TCGAA", 1, {0: {2}, 1: {2}}),
            ("TCGAA", 0, {0: {2}}),
        ],
    )
    def test_map_kmer_to_index(
        self, mapping_pipeline: MappingPipeline, kmer, n_mismatches, expected_dict
    ):
        mapping_pipeline._generate_and_write_index(1, verbose=False)
        mapping_dict = mapping_pipeline.map_kmer_to_index(kmer, n_mismatches)
        assert mapping_dict == expected_dict

    @pytest.mark.parametrize(
        "kmer,mapping_dict,expected_pos,expected_mismatches",
        [
            ("TTCGA", {0: set(), 1: {1, 2}}, [1, 2], 1),
            ("TACGA", {0: set(), 1: {2}}, [2], 1),
            ("CGAAA", {0: {3}}, [3], 0),
        ],
    )
    def test_add_kmer_to_mapping_df(
        self,
        mapping_pipeline: MappingPipeline,
        kmer,
        mapping_dict,
        expected_pos,
        expected_mismatches,
    ):
        mapping_pipeline.add_kmer_to_mapping_df(1, kmer, mapping_dict)
        result_row = mapping_pipeline.matrices["mapping"].loc[kmer]
        assert result_row["reference_sequence_position"] == expected_pos
        assert result_row["number_of_mismatches"] == expected_mismatches

    @pytest.mark.parametrize(
        "kmer,expected_group,expected_exclusivity",
        [
            ("TTCGA", "group2_enriched", "exclusive"),
            ("TACGA", "group1_enriched", "exclusive"),
            ("CGAAA", "conserved", "non-exclusive"),
        ],
    )
    def test_update_mapping_df_group(
        self,
        mapping_pipeline: MappingPipeline,
        kmer,
        expected_group,
        expected_exclusivity,
    ):
        mapping_pipeline.update_mapping_df_group(kmer)
        result_row = mapping_pipeline.matrices["mapping"].loc[kmer]
        assert result_row["group"] == expected_group
        assert result_row["exclusivity"] == expected_exclusivity

    @pytest.mark.parametrize(
        "kmer,positions,expected_gene",
        [
            ("TTCGA", [1, 2], "gene1, gene2"),
            ("TACGA", [2], "gene2"),
            ("CGAAA", [3], ""),
        ],
    )
    def test_update_mapping_df_gene(
        self,
        mapping_pipeline: MappingPipeline,
        kmer,
        positions,
        expected_gene,
    ):
        mapping_pipeline.matrices["mapping"].at[
            kmer, "reference_sequence_position"
        ] = positions
        mapping_pipeline.update_mapping_df_genes(kmer)

        result_row = mapping_pipeline.matrices["mapping"].loc[kmer]
        assert result_row["gene"] == expected_gene

    def test_fill_coverage_dataframe(
        self,
        mapping_pipeline: MappingPipeline,
        expected_coverage_df: pd.DataFrame,
    ):
        mapping_pipeline.matrices["mapping"] = pd.DataFrame(
            {
                "reference_sequence_position": [[1, 2], [2], [3]],
                "gene": ["gene1, gene2", "gene2", ""],
                "number_of_mismatches": [1, 1, 0],
                "group": [
                    "group2_enriched",
                    "group1_enriched",
                    "conserved",
                ],
                "exclusivity": [
                    "exclusive",
                    "exclusive",
                    "non-exclusive",
                ],
            },
            index=["TTCGA", "TACGA", "CGAAA"],
            columns=[
                "reference_sequence_position",
                "gene",
                "number_of_mismatches",
                "group",
                "exclusivity",
            ],
        )
        mapping_pipeline.fill_coverage_dataframe()

        result_coverage = mapping_pipeline.matrices["coverage"]
        pd.testing.assert_frame_equal(result_coverage, expected_coverage_df)

    def test_save_mapping_and_coverage_many_digits(
        self,
        mapping_pipeline: MappingPipeline,
        test_project_dir,
    ):
        mapping_pipeline.matrices["mapping"] = pd.DataFrame(
            {
                "reference_sequence_position": [[35, 154]],
                "gene": [""],
                "number_of_mismatches": [0],
                "group": [
                    "conserved",
                ],
                "exclusivity": [
                    "non-exclusive",
                ],
            },
            index=["CGAAA"],
            columns=[
                "reference_sequence_position",
                "gene",
                "number_of_mismatches",
                "group",
                "exclusivity",
            ],
        )
        mapping_pipeline.matrices["coverage"] = pd.DataFrame(
            {
                "group1_enriched": [0, 1, 1, 1, 1, 1, 0],
                "group2_enriched": [1, 2, 2, 2, 2, 1, 0],
                "conserved": [0, 0, 1, 1, 1, 1, 1],
            },
            index=[pos for pos in range(1, 8)],
        )
        os.makedirs(test_project_dir / "output", exist_ok=True)
        mapping_pipeline.save_mapping_and_coverage(1)

        assert os.path.exists(
            f"{test_project_dir}/output/mapping_5mers_max_1_mismatches.tsv"
        )
        assert os.path.exists(
            f"{test_project_dir}/output/coverage_5mers_max_1_mismatches.tsv"
        )
        saved_mapping_df = pd.read_csv(
            f"{test_project_dir}/output/mapping_5mers_max_1_mismatches.tsv",
            sep="\t",
            index_col=0,
        ).fillna("")

        assert saved_mapping_df.at["CGAAA", "reference_sequence_position"] == "35,154"

    def test_run_mapping(
        self,
        mapping_pipeline: MappingPipeline,
        test_project_dir: Path,
        expected_mapping_df: pd.DataFrame,
        expected_coverage_df: pd.DataFrame,
    ):
        os.makedirs(test_project_dir / "output", exist_ok=True)

        mapping_pipeline.run_mapping(1, verbose=False)

        result_mapping_df = pd.read_csv(
            f"{test_project_dir}/output/mapping_5mers_max_1_mismatches.tsv",
            sep="\t",
            index_col=0,
        ).fillna("")

        result_coverage_df = pd.read_csv(
            f"{test_project_dir}/output/coverage_5mers_max_1_mismatches.tsv",
            sep="\t",
            index_col=0,
        ).fillna("")

        pd.testing.assert_frame_equal(result_mapping_df, expected_mapping_df)
        pd.testing.assert_frame_equal(result_coverage_df, expected_coverage_df)
