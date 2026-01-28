import pandas as pd
import os
import yaml
import pathlib
import pickle
from Bio import SeqIO


class BasePipeline:
    def _load_and_check_config(self, project_dir: str) -> None:
        try:
            with open(f"{project_dir}/config.yaml", "r") as config_file:
                config = yaml.safe_load(config_file)
            if len(config["sample_groups"]) != len(config["input_files"]):
                err_msg = "Number of sample groups and input fastas must match!"
                raise ValueError()
            if len(config["header_format"]) == 0:
                self.header_format = [
                    "plain" for _ in range(len(config["sample_groups"]))
                ]
            else:
                if len(config["header_format"]) != len(config["sample_groups"]):
                    err_msg = "Number of header formats must be 0 or match number of sample groups!"
                    raise ValueError()
                if any(
                    [
                        header != "gisaid" and header != "ncbi" and header != "plain"
                        for header in config["header_format"]
                    ]
                ):
                    err_msg = "Header format must be 'gisaid', 'ncbi' or 'plain'!"
                    raise ValueError()
                self.header_format = config["header_format"]
            self.sample_groups = config["sample_groups"]
            self.input_files = config["input_files"]
            self.mode = config["mode"]
            if self.mode == "ovr":
                self.goi = config["goi"]
                self.control_groups = self.sample_groups.copy()
                self.control_groups.remove(self.goi)
            elif self.mode == "ava":
                self.goi = None
                self.control_groups = None
            else:
                err_msg = "Unrecognized analysis mode, should by ava or ovr. Check your config file!"
                raise ValueError()
            self.ref_seq_name = ".".join(config["ref_seq"].split(".")[:-1])
            try:
                self.ref_seq = str(
                    SeqIO.read(
                        f"{project_dir}/input/{config['ref_seq']}", format="fasta"
                    ).seq
                )
            except ValueError:
                err_msg = (
                    "Missing reference sequence file or the file has incorrect format!"
                )
                raise ValueError()
            try:
                self.ref_genes = config["ref_genes"]
            except KeyError:
                self.ref_genes = {}

        except FileNotFoundError:
            raise FileNotFoundError(
                f"No config.yaml file found in directory {project_dir}!"
            )
        except ValueError:
            raise ValueError(err_msg)
        except KeyError:
            raise KeyError("Config file is missing required fields!")

    def __init__(self, project_dir: str, k: int = 5):
        if os.path.isdir(project_dir) == True:
            self.project_dir = project_dir
        else:
            raise NotADirectoryError("Project directory does not exist!")

        if type(k) != int:
            raise TypeError("k must by an integer!")
        else:
            self.k = k

        self.sample_groups = []
        self.input_files = []
        self.header_format = []
        self.mode = ""
        self.goi = ""
        self.control_groups = []
        self.kmer_set = set()
        self.sample_list = []
        self.sample_group_dict = {}
        self.ref_seq_name = ""
        self.ref_seq = ""
        self.ref_genes = {}

        self._load_and_check_config(project_dir)

    def load_kmer_set(self, k) -> list:
        kmer_set = set()
        kmer_set_in_path = list(
            pathlib.Path(f"{self.project_dir}/input/indices/").glob(f"{k}mer_*_set.pkl")
        )
        if len(kmer_set_in_path) != len(self.sample_groups):
            raise FileNotFoundError(
                f"Some or all {k}-mers are missing from the project directory! Please run PORTEK find_k!"
            )
        for filename in kmer_set_in_path:
            with open(filename, mode="rb") as in_file:
                partial_set = pickle.load(in_file)
            kmer_set.update(partial_set)
        kmer_set = list(kmer_set)
        self.kmer_set = kmer_set
        return kmer_set

    def load_sample_list(self) -> list:
        sample_list = []
        sample_list_in_path = list(
            pathlib.Path(f"{self.project_dir}/input/indices").glob("*sample_list.pkl")
        )
        if len(sample_list_in_path) != len(self.sample_groups):
            raise FileNotFoundError(
                "Some or all samples are missing from the project directory! Please run PORTEK find_k!"
            )
        for filename in sample_list_in_path:
            with open(filename, mode="rb") as in_file:
                partial_list = pickle.load(in_file)
            group = filename.stem.split("_")[0]
            partial_list = [f"{group}_{sample_name}" for sample_name in partial_list]
            sample_list.extend(partial_list)
        sample_group_dict = {
            f"{group}": [
                sample for sample in sample_list if sample.split("_")[0] == f"{group}"
            ]
            for group in self.sample_groups
        }

        self.sample_list = sample_list
        self.sample_group_dict = sample_group_dict
        return sample_list

    def check_min_max_k(self, mink: int, maxk: int) -> None:
        if type(mink) != int or mink < 5 or mink % 2 == 0:
            raise TypeError("Minimum k must by an odd integer not smaller than 5!")
        else:
            self.mink = mink
        if type(maxk) != int or maxk < 5 or maxk % 2 == 0:
            raise TypeError("Maximum k must by an odd integer not smaller than 5!")
        else:
            self.maxk = maxk
        if self.maxk < self.mink:
            raise ValueError("Minimum k must be no greater than maximum k!")


def encode_seq(nuc_seq: str) -> list[str]:
    nuc_seq_upper = nuc_seq.upper()
    encoding = {"A": "00", "C": "01", "G": "10", "T": "11"}
    bit_seq = [encoding.get(nuc, "X") for nuc in nuc_seq_upper]
    return bit_seq


def decode_kmer(id: int, k) -> str:
    decoding = {"00": "A", "01": "C", "10": "G", "11": "T"}
    kmer_bin_string = bin(id)[2:].rjust(2 * k, "0")
    kmer_bin_string = [
        kmer_bin_string[i : i + 2] for i in range(0, len(kmer_bin_string), 2)
    ]
    kmer_seq = "".join([decoding[bits] for bits in kmer_bin_string])
    return kmer_seq


def assign_kmer_group_ava(
    row: pd.Series, p_cols: list, avg_cols: list, freq_cols: list, err_cols: list
):
    max_group = row[avg_cols].idxmax()
    max_group = max_group.split("_")[0]
    rel_p_cols = [col for col in p_cols if max_group in col]
    if all(row[freq_cols] > 0.9) and all((abs(row[err_cols]) < 0.1)):
        return "conserved"
    elif all(row[rel_p_cols] < 0.01):
        return f"{max_group}_enriched"
    else:
        return "not_significant"


def assign_kmer_group_ovr(
    row: pd.Series, goi: str, p_cols: list, err_cols: list, freq_cols: list
):
    if all(row[freq_cols] > 0.9) and all((abs(row[err_cols]) < 0.1)):
        return "conserved"
    elif all(row[p_cols] < 0.01) and all(row[err_cols] > 0):
        return f"{goi}_enriched"
    elif all(row[p_cols] < 0.01) and all(row[err_cols] < 0):
        return "control_enriched"
    elif any(row[p_cols] < 0.01):
        return "group_dependent"

    else:
        return "not_significant"


def check_exclusivity(row: pd.Series, avg_cols: list) -> str:
    if len([col for col in avg_cols if row[col] > 0]) == 1:
        return "exclusive"
    else:
        return "non-exclusive"
