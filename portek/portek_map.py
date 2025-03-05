import itertools
import math
import operator
import os
import pathlib
import pickle
import shutil
import subprocess

import numpy as np
import pandas as pd
import pysam
import regex
import yaml
from Bio import SeqIO
from scipy.ndimage import histogram
from scipy.signal import find_peaks

import portek


class MappingPipeline:

    def __init__(self, project_dir: str, k):
        if os.path.isdir(project_dir) == True:
            self.project_dir = project_dir
        else:
            raise NotADirectoryError("Project directory does not exist!")

        if type(k) != int:
            raise TypeError("k must by an integer!")
        else:
            self.k = k

        try:
            with open(f"{project_dir}/config.yaml", "r") as config_file:
                config = yaml.safe_load(config_file)
            self.sample_groups = config["sample_groups"]
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
                err_msg = "No or wrong reference sequence file!"
                raise ValueError()
            if "ref_genes" in config.keys():
                self.ref_genes = config["ref_genes"]
            else:
                self.ref_genes = None

            self.avg_cols = [f"{group}_avg" for group in self.sample_groups]

        except FileNotFoundError:
            raise FileNotFoundError(
                f"No config.yaml file found in directory {project_dir} or the file has missing/wrong configuration!"
            )
        except ValueError:
            raise ValueError(err_msg)

        self.matrices = {}
        try:
            self.matrices["enriched"] = pd.read_csv(
                f"{project_dir}/output/enriched_{self.k}mers_stats.csv", index_col=0
            )

        except:
            raise FileNotFoundError(
                f"No enriched {self.k}-mers table found in {project_dir}output/ ! Please run PORT-EK enriched first!"
            )

        self.mutations = None
        self.sample_list = None
        self.sample_group_dict = None
        self.group_avg_pos = None

    def get_samples(self, verbose: bool = False):
        sample_list_in_path = pathlib.Path(f"{self.project_dir}/input/indices").glob(
            "*sample_list.pkl"
        )
        sample_list = []
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
        if len(sample_list) == 0:
            raise FileNotFoundError("No sample lists found!")
        self.sample_list = sample_list
        self.sample_group_dict = sample_group_dict

    def _check_bowtie2_path(self):
        return shutil.which("bowtie2")

    def _check_index_built(self):
        index_files = list(
            pathlib.Path(f"{self.project_dir}/temp/ref_index/").glob(
                f"{self.ref_seq_name}.*"
            )
        )
        if len(index_files) == 0:
            return False
        else:
            return True

    def _bowtie_build_index(self, verbose: bool = False):
        if os.path.exists(f"{self.project_dir}/temp/ref_index/") == False:
            os.makedirs(f"{self.project_dir}/temp/ref_index")
        build_cmd = [
            f"bowtie2-build",
            "-f",
            f"{self.project_dir}/input/{self.ref_seq_name}.fasta",
            f"{self.project_dir}/temp/ref_index/{self.ref_seq_name}",
        ]
        result = subprocess.run(build_cmd, capture_output=True, text=True)
        if verbose == True:
            print(build_cmd)
        if result.returncode != 0:
            raise Exception(result.stderr)
        else:
            if verbose == True:
                print(result.stdout)

    def _bowtie_map(self, verbose: bool = False):
        seed_length = int(math.ceil(self.k / 2))
        map_cmd = [
            f"bowtie2",
            "--norc",
            "-a",
            "-L",
            f"{seed_length}",
            "-x",
            f"{self.project_dir}/temp/ref_index/{self.ref_seq_name}",
            "-f",
            f"{self.project_dir}/temp/enriched_{self.k}mers.fasta",
            "-S",
            f"{self.project_dir}/temp/enriched_{self.k}mers.sam",
        ]
        if verbose == True:
            print(" ".join(map_cmd))
        result = subprocess.run(map_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise Exception(result.stderr)
        else:
            if verbose == True:
                print(result.stdout)

    def run_mapping(self, verbose: bool = False):
        if self._check_bowtie2_path() == None:
            raise FileNotFoundError(
                f"bowtie2 not found! Please install bowtie and add it to your PATH!"
            )
        if self._check_index_built() == False:
            self._bowtie_build_index(verbose=verbose)
        self._bowtie_map(verbose=verbose)

    def _parse_CIGAR(self, CIGAR_string: str) -> dict:
        n_repeats = regex.findall(r"\d+", CIGAR_string)
        matches = [char for char in CIGAR_string if char.isalpha()]
        CIGAR_list = []
        for n, m in zip(n_repeats, matches):
            CIGAR_list.extend(int(n) * [m])
        return CIGAR_list

    def _read_sam_to_df(self) -> pd.DataFrame:
        reads = pysam.AlignmentFile(
            f"{self.project_dir}/temp/enriched_{self.k}mers.sam", mode="r"
        )

        read_dict = {
            "kmer": [],
            "flag": [],
            "ref_pos": [],
            "CIGAR": [],
            "n_mismatch": [],
            "score": [],
        }
        for read in reads:
            read_dict["kmer"].append(read.query_name)
            read_dict["flag"].append(read.flag)
            if read.reference_start == -1:
                read_dict["ref_pos"].append(0)
            else:
                read_dict["ref_pos"].append(read.reference_start)
            if read.cigarstring == None:
                read_dict["CIGAR"].append("")
            else:
                read_dict["CIGAR"].append(read.cigarstring)
            if read.has_tag("NM"):
                read_dict["n_mismatch"].append(read.get_tag("NM", False))
            else:
                read_dict["n_mismatch"].append(0)
            if read.has_tag("AS"):
                read_dict["score"].append(read.get_tag("AS", False))
            else:
                read_dict["score"].append(-6 * self.k)

        mappings_df = pd.DataFrame(read_dict)
        mappings_df["CIGAR"] = mappings_df["CIGAR"].apply(self._parse_CIGAR)
        mappings_df.loc[:, ["flag", "ref_pos", "n_mismatch"]] = mappings_df.loc[
            :, ["flag", "ref_pos", "n_mismatch"]
        ].astype(int)
        mappings_df["group"] = self.matrices["enriched"]["group"]
        mappings_df["mutations"] = "WT"
        mappings_df["mapping_ok"] = 1
        mappings_df["ref_pos"] = mappings_df.apply(
            lambda row: row["ref_pos"] + 1 if row["flag"] != 4 else row["ref_pos"], axis=1
        )
        mappings_df["group"] = mappings_df["kmer"].apply(
            lambda kmer: self.matrices["enriched"].loc[kmer, "group"]
        )
        return mappings_df

    def _align_seqs(self, ref_seq, kmer, map_pos, cigar):
        ref_start = map_pos
        aln_len = len(cigar)
        ref_end = ref_start + aln_len
        q_seq = [nuc for nuc in kmer]
        t_seq = [nuc for nuc in ref_seq[ref_start:ref_end]]
        ref_pos = []
        curr_pos = ref_start - 1
        for i, change in enumerate(cigar):
            if change == "D":
                curr_pos += 1
                q_seq.insert(i, "-")
                ref_pos.append(curr_pos)
            elif change == "I":
                t_seq.insert(i, "-")
                ref_pos.append(curr_pos)
            else:
                curr_pos += 1
                ref_pos.append(curr_pos)

        q_seq = q_seq[:aln_len]
        t_seq = t_seq[:aln_len]
        return q_seq, t_seq, ref_pos

    def _join_indels(self, mutations: list):
        subs = []
        ins_pos = []
        del_pos = []

        for mut in mutations:
            if mut[1] == "ins":
                ins_pos.append(mut[0])
            elif mut[1] == "del":
                del_pos.append(mut[0])
            else:
                subs.append(mut)

        ins_pos = set(ins_pos)
        inss = []
        dels = []

        for start_pos in ins_pos:
            muts = [
                mut[2] for mut in mutations if mut[0] == start_pos and mut[1] == "ins"
            ]
            inss.append((start_pos, "ins", "".join(muts)))

        for k, g in itertools.groupby(
            enumerate(sorted(del_pos)), lambda x: x[0] - x[1]
        ):
            del_group_pos = list(map(operator.itemgetter(1), g))
            dels.append((del_group_pos[0], "del", del_group_pos[-1]))
        grouped_muts = subs + inss + dels
        grouped_muts.sort()
        return grouped_muts

    def _find_variants(self, ref_seq, kmer, map_pos, cigar) -> list:
        map_pos = map_pos-1
        mutations = []
        q_seq, t_seq, ref_pos = self._align_seqs(ref_seq, kmer, map_pos, cigar)
        if len(q_seq) != len(t_seq) != len(ref_pos):
            raise ValueError(
                f"Improper alingment of k-mer {q_seq}, reference {t_seq}, and refrence position {ref_pos}"
            )
        for i in range(len(q_seq)):
            if q_seq[i] != t_seq[i]:
                if q_seq[i] == "-":
                    mutations.append((ref_pos[i] + 1, "del", ref_pos[i] + 1))
                elif t_seq[i] == "-":
                    mutations.append((ref_pos[i] + 1, "ins", q_seq[i]))
                else:
                    mutations.append((ref_pos[i] + 1, t_seq[i], q_seq[i]))

        mutations = self._join_indels(mutations)
        return mutations

    def _mutation_tuple_to_text(self, mutation_as_tuple: tuple) -> str:
        if mutation_as_tuple[1] == "del":
            mutation_as_text = f"{mutation_as_tuple[0]}_{mutation_as_tuple[2]}del"
        elif mutation_as_tuple[1] == "ins":
            mutation_as_text = f"{mutation_as_tuple[0]}_{mutation_as_tuple[0]+1}ins{mutation_as_tuple[2]}"
        else:
            mutation_as_text = (
                f"{mutation_as_tuple[0]}{mutation_as_tuple[1]}>{mutation_as_tuple[2]}"
            )
        return mutation_as_text

    def _load_kmer_pos(
        self,
        input_filename: pathlib.Path,
        kmers: np.ndarray,
        kmer_distros: dict,
        verbose: bool = False,
    ) -> None:
        with open(input_filename, mode="rb") as in_file:
            group = f"{input_filename.stem.split('_')[1]}_enriched"
            temp_dict = pickle.load(in_file)
        temp_dict = {
            portek.decode_kmer(id, self.k): positions
            for id, positions in temp_dict.items()
        }
        kmer_distros[group] = {
            kmer: positions for kmer, positions in temp_dict.items() if kmer in kmers
        }

    def _get_total_distros(self, kmer_distros: dict) -> None:
        total_distros = {}
        for distros in kmer_distros.values():
            for kmer, distro in distros.items():
                if kmer in total_distros.keys():
                    total_distros[kmer] = total_distros[kmer] + distro
                else:
                    total_distros[kmer] = distro
        kmer_distros["conserved"] = total_distros

    def _get_peaks_from_distro(self, kmer_distros: dict) -> dict:
        distro_peaks = {group: {} for group in kmer_distros.keys()}
        for group in distro_peaks.keys():
            for kmer, distro in kmer_distros[group].items():
                max_pos = max(max(distro), len(self.ref_seq))
                bins = max_pos + 1
                histogram_over_genome = histogram(distro, 0, max_pos, bins)
                peaks = find_peaks(histogram_over_genome, distance=bins / 10)
                distro_peaks[group][kmer] = peaks[0]
        return distro_peaks

    def _get_kmer_peaks(self, mappings_df: pd.DataFrame, verbose: bool = False):
        print(f"\nLoading enriched {self.k}-mer position distributions.")
        in_path = pathlib.Path(f"{self.project_dir}/input/indices/").glob(
            f"{self.k}mer_*_pos_dict.pkl"
        )
        if verbose == True:
            files_counter = 1
            tot_files = len(self.sample_groups)
        kmer_distros = {}
        kmers = mappings_df["kmer"].unique()
        for filename in in_path:
            self._load_kmer_pos(filename, kmers, kmer_distros, verbose)
            if verbose == True:
                print(
                    f"Loaded {self.k}-mer distributions from {files_counter} of {tot_files} groups.",
                    end="\r",
                    flush=True,
                )
                files_counter += 1
        self._get_total_distros(kmer_distros)
        distro_peaks = self._get_peaks_from_distro(kmer_distros)
        return distro_peaks

    def _get_real_pos(
        self, group: str, kmer: str, ref_pos: int, actual_positions: dict
    ) -> int:
        peaks = actual_positions[group][kmer]
        diffs = abs(peaks - ref_pos)
        return peaks[np.argmin(diffs)]

    def _get_n_peaks(
        self, group: str, kmer: str, ref_pos: int, actual_positions: dict
    ) -> int:
        peaks = actual_positions[group][kmer]
        return len(peaks)

    def _resolve_multiple_mappings(self, mappings_df: pd.DataFrame) -> pd.Index:
        multimapped_kmers = np.unique(
            mappings_df[mappings_df["flag"] == 256]["kmer"].to_numpy()
        )
        bad_mappings = pd.Index([])
        for kmer in multimapped_kmers:
            kmer_sub_df = mappings_df.loc[mappings_df["kmer"] == kmer].copy()
            kmer_sub_df["pos_diff"] = abs(
                kmer_sub_df["ref_pos"] - kmer_sub_df["real_pos"]
            )
            kmer_sub_df.sort_values(
                ["score", "pos_diff"], ascending=[False, True], inplace=True
            )
            n_peaks = kmer_sub_df["n_peaks"].values[0]
            bad_mappings = bad_mappings.union(kmer_sub_df.iloc[n_peaks:].index)
        return bad_mappings

    def _filter_mappings(self, mappings_df: pd.DataFrame):
        bad_mappings = self._resolve_multiple_mappings(mappings_df)
        mappings_df.loc[bad_mappings, "mapping_ok"] = 0
        mappings_df = mappings_df[mappings_df["mapping_ok"] == 1]
        return mappings_df

    def _format_mappings_df(self, mappings_df: pd.DataFrame):
        mappings_df.loc[mappings_df["flag"] == 4, "mutations"] = "-"
        mappings_df.loc[mappings_df["flag"] == 4, "n_mismatch"] = self.k
        formatted_df = mappings_df.drop(
            ["flag", "CIGAR", "score", "mapping_ok", "real_pos", "n_peaks"],
            axis=1
        ).sort_values("ref_pos")
        return formatted_df

    def _count_mappings(self, mappings_df: pd.DataFrame):
        num_kmers = len(self.matrices["enriched"])
        validly_aln_kmers = mappings_df[mappings_df["flag"] != 4]["kmer"].to_numpy()
        num_aligned = len(np.unique(validly_aln_kmers))
        num_multiple = len(validly_aln_kmers) - num_aligned
        num_unaligned = num_kmers - num_aligned
        print(
            f"\n{num_aligned} out of {num_kmers} {self.k}-mers were aligned to the reference."
        )
        print(f"Of those, {num_multiple} were aligned to more than one position.")
        print(f"{num_unaligned} {self.k}-mers couldn't be aligned.")

    def analyze_mapping(self, verbose: bool = False):
        mappings_df = self._read_sam_to_df()
        actual_positions = self._get_kmer_peaks(mappings_df, verbose)
        mappings_df["real_pos"] = mappings_df.apply(
            lambda row: self._get_real_pos(
                row["group"], row["kmer"], row["ref_pos"], actual_positions
            ),
            axis=1,
        )
        mappings_df["n_peaks"] = mappings_df.apply(
            lambda row: self._get_n_peaks(
                row["group"], row["kmer"], row["ref_pos"], actual_positions
            ),
            axis=1,
        )
        filtered_df = self._filter_mappings(mappings_df)
        for row in filtered_df.itertuples():
            if row.n_mismatch > 0:
                mutations_as_tuples = self._find_variants(
                    self.ref_seq, row.kmer, row.ref_pos, row.CIGAR
                )
                filtered_df.loc[row.Index, "mutations"] = "; ".join(
                    [self._mutation_tuple_to_text(mut) for mut in mutations_as_tuples]
                )

        if verbose == True:
            self._count_mappings(filtered_df)

        formatted_df = self._format_mappings_df(filtered_df)
        self.matrices["mappings"] = formatted_df


    def save_mappings_df(self):
        print(f"\nSaving {self.k}-mer mappings.")
        df_to_save = self.matrices["mappings"].copy()
        df_to_save.index.name = "id"
        df_to_save.to_csv(
            f"{self.project_dir}/output/enriched_{self.k}mers_mappings.csv"
        )
