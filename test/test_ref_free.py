import pytest
import os
import pandas as pd
from unittest.mock import patch, mock_open, MagicMock
from portek.portek_map import RefFreePipeline


class TestRefFreePipeline:

    @patch("os.path.isdir", return_value=False)
    def test_constructor_invalid_directory(self, mock_isdir):
        with pytest.raises(NotADirectoryError):
            RefFreePipeline("invalid_dir", 15)

    @patch("os.path.isdir", return_value=True)
    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data="sample_groups: ['group1', 'group2']\nmode: 'ava'",
    )
    @patch("pandas.read_csv")
    def test_constructor_valid_directory(self, mock_read_csv, mock_open, mock_isdir):
        mock_read_csv.return_value = pd.DataFrame(
            {"group": ["group1", "group2"], "index": ["AAA", "CCC"]}
        ).set_index("index")

        pipeline = RefFreePipeline("valid_dir", 15)
        assert pipeline.project_dir == "valid_dir"
        assert pipeline.k == 15
        assert pipeline.sample_groups == ["group1", "group2"]
        assert pipeline.mode == "ava"

    @patch("os.path.isdir", return_value=True)
    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data="sample_groups: ['group1', 'group2']\nmode: 'ovr'\ngoi: 'group1'",
    )
    @patch("pandas.read_csv")
    def test_constructor_ovr_mode(self, mock_read_csv, mock_open, mock_isdir):
        mock_read_csv.return_value = pd.DataFrame(
            {"group": ["group1", "group2"], "index": ["AAA", "CCC"]}
        ).set_index("index")

        pipeline = RefFreePipeline("valid_dir", 15)
        assert pipeline.goi == "group1"
        assert pipeline.control_groups == ["group2"]

    @patch("os.path.isdir", return_value=True)
    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data="sample_groups: ['group1', 'group2']\nmode: 'invalid_mode'",
    )
    def test_constructor_invalid_mode(self, mock_open, mock_isdir):
        with pytest.raises(ValueError):
            RefFreePipeline("valid_dir", 15)

    @patch("os.path.isdir", return_value=True)
    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data="sample_groups: ['group1', 'group2']\nmode: 'ava'",
    )
    def test_constructor_missing_config_file(self, mock_open, mock_isdir):
        mock_open.side_effect = FileNotFoundError
        with pytest.raises(FileNotFoundError):
            RefFreePipeline("valid_dir", 15)

    @patch("os.path.isdir", return_value=True)
    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data="sample_groups: ['group1', 'group2']\nmode: 'ava'",
    )
    @patch("pandas.read_csv", side_effect=FileNotFoundError)
    def test_constructor_missing_enriched_file(
        self, mock_read_csv, mock_open, mock_isdir
    ):
        with pytest.raises(FileNotFoundError):
            RefFreePipeline("valid_dir", 15)

    @patch("portek.portek_map.pd.DataFrame.to_csv")
    @patch("os.path.isdir", return_value=True)
    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data="sample_groups: ['group1', 'group2']\nmode: 'ava'",
    )
    def test_save_group_distros_saves_csv(self, mock_open, mock_isdir, mock_to_csv):
        # Setup
        pipeline = RefFreePipeline("valid_dir", 15)
        pipeline.group_avg_pos = {"pos1": [1, 2], "pos2": [3, 4]}
        pipeline.matrices = {"enriched": pd.DataFrame({"group": ["A", "B"]})}
        pipeline.project_dir = "/valid/path"

        # Execute
        pipeline.save_group_distros()

        # Verify
        mock_to_csv.assert_called_once_with("/valid/path/temp/group_avg_pos.csv")

    @patch("os.path.isdir", return_value=True)
    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data="sample_groups: ['group1', 'group2']\nmode: 'ava'",
    )
    def test_save_group_distros_missing_group_column(
        self,
        mock_open,
        mock_isdir,
    ):
        # Setup
        pipeline = RefFreePipeline("valid_dir", 15)
        pipeline.group_avg_pos = {"pos1": [1, 2], "pos2": [3, 4]}
        pipeline.matrices = {"enriched": pd.DataFrame()}  # Missing 'group' column
        pipeline.project_dir = "/valid/path"

        # Execute & Verify
        with pytest.raises(KeyError):
            pipeline.save_group_distros()
