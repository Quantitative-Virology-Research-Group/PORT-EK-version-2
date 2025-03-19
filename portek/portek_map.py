import itertools
import math
import multiprocessing
import operator
import os
import pathlib
import pickle
import shutil
import subprocess
import json

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import regex
import seaborn as sns
import yaml
from Bio import SeqIO
from scipy.ndimage import histogram
from scipy.signal import find_peaks
from scipy.stats import linregress

import portek
from portek.portek_utils import BasePipeline


class MappingPipeline(BasePipeline):

    def __init__(self, project_dir: str, k: int):
        super().__init__(project_dir, k)
        self.kmer_index = None
        self.matrices = {}

        try:
            self.matrices["enriched"] = pd.read_csv(
                f"{project_dir}/output/enriched_{self.k}mers_stats.csv", index_col=0
            )
        except:
            raise FileNotFoundError(
                f"No enriched {self.k}-mers table found in {project_dir}output/ ! Please run PORT-EK find_enriched first!"
            )

    def _check_index(self, dist: int) -> bool:
        return pathlib.Path(
            f"{self.project_dir}/temp/ref_index/{self.ref_seq_name}_index_{self.k}_{dist}.pkl"
        ).exists()

    def index_ref_seq(self, dist: int, verbose: bool = False) -> None:
        if self._check_index(dist) == False:
            if verbose == True:
                print(
                    f"No {self.k}-mer index with maximum distance {dist} exists for {self.ref_seq_name}, building index."
                )
            bit_ref_seq = portek.encode_seq(self.ref_seq)
            kmer_index = {dist: {} for dist in range(dist + 1)}
            for i in range(0, len(bit_ref_seq) - self.k + 1):
                bit_kmer = bit_ref_seq[i : i + self.k]
                if "X" not in bit_kmer:
                    int_kmer0 = int("".join(bit_kmer), base=2)
                    if int_kmer0 in kmer_index[0].keys():
                        kmer_index[0][int_kmer0].append(i + 1)
                    else:
                        kmer_index[0][int_kmer0] = [i + 1]

                    for d in range(1, dist + 1):
                        all_bit_kmers_del_d = set(
                            itertools.combinations(bit_kmer, len(bit_kmer) - d)
                        )
                        all_int_kmers_del_d = [
                            int("".join(bit_kmer_del_d), base=2)
                            for bit_kmer_del_d in all_bit_kmers_del_d
                        ]
                        for int_kmer_del_d in all_int_kmers_del_d:
                            if int_kmer_del_d in kmer_index[d].keys():
                                kmer_index[d][int_kmer_del_d].append(i + 1)
                            else:
                                kmer_index[d][int_kmer_del_d] = [i + 1]

            with open(
                f"{self.project_dir}/temp/ref_index/{self.ref_seq_name}_index_{self.k}_{dist}.pkl",
                mode="wb",
            ) as out_file:
                pickle.dump(kmer_index, out_file)
            if verbose == True:
                print(
                    f"Finished building {self.k}-mer index with maximum distance {dist} exists for {self.ref_seq_name}"
                )
                for d in range(dist + 1):
                    print(f"Extracted {len(kmer_index[d])} k-mers with distance {d}.")
            self.kmer_index = kmer_index
        else:
            if verbose == True:
                print(
                    f"{self.k}-mer index with maximum distance {dist} for {self.ref_seq_name} already exists."
                )
            with open(
                f"{self.project_dir}/temp/ref_index/{self.ref_seq_name}_index_{self.k}_{dist}.pkl",
                mode="rb",
            ) as in_file:
                kmer_index = pickle.load(in_file)
            self.kmer_index = kmer_index

    def run_mapping(self, dist: int, verbose: bool = False) -> dict:
        mapping_dict = {kmer:{} for kmer in self.matrices["enriched"].index}
        unmapped_counter = {d:0 for d in range(dist+1)}
        for kmer in mapping_dict.keys():
            for d in range(dist+1):
                mapping_dict[kmer][d]=[]
                bit_kmer = portek.encode_seq(kmer)
                bit_kmers_to_check = set(
                            itertools.combinations(bit_kmer, len(bit_kmer) - d)
                        )
                int_kmers_to_check = [
                            int("".join(bit_kmer_del_d), base=2)
                            for bit_kmer_del_d in bit_kmers_to_check
                        ]
                for int_kmer in int_kmers_to_check:
                    mapping_dict[kmer][d].extend(self.kmer_index[d].get(int_kmer,[]))
                if len(mapping_dict[kmer][d]) == 0:
                    unmapped_counter[d] += 1
        print(len(mapping_dict.keys()))
        print(unmapped_counter)
        return mapping_dict
    
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
        map_pos = map_pos - 1
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
        self, input_filename: pathlib.Path, kmers: np.ndarray, verbose: bool = False
    ) -> None:
        if verbose == True:
            print(f"Loading k-mer positions from file {input_filename.stem}.")
        with open(input_filename, mode="rb") as in_file:
            group = f"{input_filename.stem.split('_')[1]}_enriched"
            temp_dict = pickle.load(in_file)
        temp_dict = {
            portek.decode_kmer(id, self.k): positions
            for id, positions in temp_dict.items()
        }
        distro = {kmer: temp_dict[kmer] for kmer in kmers if kmer in temp_dict.keys()}
        if verbose == True:
            print(f"Done loading k-mer positions from file {input_filename.stem}.")
        return group, distro

    def _get_group_distros(self, pos_results: list) -> dict:
        kmer_distros = {}
        if self.mode == "ava":
            for result in pos_results:
                kmer_distros[result[0]] = result[1]
        else:
            kmer_distros["control_enriched"] = {}
            for result in pos_results:
                if result[0] == f"{self.goi}_enriched":
                    kmer_distros[result[0]] = result[1]
                else:
                    kmer_distros["control_enriched"] = {
                        kmer: kmer_distros["control_enriched"].get(kmer, [])
                        + result[1].get(kmer, [])
                        for kmer in set(
                            list(kmer_distros["control_enriched"].keys())
                            + list(result[1].keys())
                        )
                    }
        return kmer_distros

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

    def _get_kmer_peaks(
        self, mappings_df: pd.DataFrame, n_jobs: int = 4, verbose: bool = False
    ):
        print(f"\nLoading enriched {self.k}-mer position distributions.")
        in_path = list(
            pathlib.Path(f"{self.project_dir}/input/indices/").glob(
                f"{self.k}mer_*_pos_dict.pkl"
            )
        )

        kmers = mappings_df["kmer"].unique()
        pool_args = []
        for filename in in_path:
            pool_args.append((filename, kmers, verbose))

        with multiprocessing.get_context("forkserver").Pool(n_jobs) as pool:
            results = pool.starmap(self._load_kmer_pos, pool_args, chunksize=1)

        kmer_distros = self._get_group_distros(results)
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

    def _filter_mappings(
        self, mappings_df: pd.DataFrame, verbose: bool = False
    ) -> pd.DataFrame:
        bad_mappings = self._resolve_multiple_mappings(mappings_df)
        mappings_df.loc[bad_mappings, "mapping_ok"] = 0
        mappings_df = mappings_df[mappings_df["mapping_ok"] == 1]
        if verbose == True:
            print(f"\nRejected {len(bad_mappings)} dubious alignments.")
        return mappings_df

    def _predict_unmapped(self, mappings_df: pd.DataFrame) -> pd.DataFrame:
        mappings_with_pred_df = mappings_df.copy()
        mappings_with_pred_df["pred_pos"] = 0
        mappings_with_pred_df["pred_err"] = 0.0
        mappings_with_pred_df["pred_r2"] = 0.0
        for group in mappings_with_pred_df["group"].unique():
            group_mappings = mappings_with_pred_df.loc[
                mappings_with_pred_df["group"] == group
            ].index
            mapped = mappings_with_pred_df.loc[
                (mappings_with_pred_df["ref_pos"] != 0)
                & (mappings_with_pred_df["group"] == group)
            ].index
            if len(mapped) < 30:
                print(
                    f"Not enough mapped k-mers for position prediction in group {group}, skipping."
                )
                continue
            real_mapped_pos = mappings_with_pred_df.loc[mapped, "real_pos"]
            ref_mapped_pos = mappings_with_pred_df.loc[mapped, "ref_pos"]
            regress = linregress(real_mapped_pos, ref_mapped_pos)
            mappings_with_pred_df.loc[group_mappings, "pred_pos"] = round(
                regress.slope * mappings_with_pred_df.loc[group_mappings, "real_pos"]
                + regress.intercept
            ).astype(int)
            pred_mapped_err = (
                mappings_with_pred_df.loc[mapped, "pred_pos"] - ref_mapped_pos
            )
            rmse = np.sqrt(sum(pred_mapped_err**2) / len(pred_mapped_err))
            mappings_with_pred_df.loc[group_mappings, "pred_err"] = round(rmse, 2)
            mappings_with_pred_df.loc[group_mappings, "pred_r2"] = round(
                regress.rvalue, 2
            )
        return mappings_with_pred_df

    def _format_mappings_df(self, mappings_df: pd.DataFrame) -> pd.DataFrame:
        mappings_df.loc[mappings_df["flag"] == 4, "mutations"] = "-"
        mappings_df.loc[mappings_df["flag"] == 4, "n_mismatch"] = self.k
        formatted_df = mappings_df.drop(
            ["flag", "CIGAR", "score", "mapping_ok", "n_peaks", "real_pos"], axis=1
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

    def analyze_mapping(self, n_jobs: int = 4, verbose: bool = False):
        mappings_df = self._read_sam_to_df()
        actual_positions = self._get_kmer_peaks(mappings_df, n_jobs, verbose)
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
        filtered_df = self._filter_mappings(mappings_df, verbose)
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
        mappings_with_pred_df = self._predict_unmapped(filtered_df)
        formatted_df = self._format_mappings_df(mappings_with_pred_df)
        self.matrices["mappings"] = formatted_df

    def _set_histogram_ax_properties(
        self,
        subplot: matplotlib.axes.Axes,
        max_pos: int,
        group: str,
        i: int,
        fig_cols: int,
    ) -> None:
        subplot.set_title(group)
        subplot.set_xlim(0, max_pos)
        subplot.set_xlabel("Reference genome position")
        if i % fig_cols == 0:
            subplot.set_ylabel("K-mer counts")

    def plot_kmer_histograms(self) -> None:
        groups = self.matrices["mappings"]["group"].unique()
        n_figs = len(groups)
        fig_cols = min(n_figs, 3)
        fig_rows = int(math.ceil(n_figs / fig_cols))
        fig, axes = plt.subplots(
            fig_rows, fig_cols, figsize=(fig_cols * 6, fig_rows * 6)
        )
        plt.subplots_adjust(hspace=0.25, wspace=0.25)
        max_pos = len(self.ref_seq)
        bins = 100
        for i, group in enumerate(groups):
            data = self.matrices["mappings"].loc[
                (self.matrices["mappings"]["group"] == group)
                & (self.matrices["mappings"]["ref_pos"] != 0),
                "ref_pos",
            ]
            if len(axes.shape) == 1:
                subplot = axes[i]
            else:
                subplot = axes[i // fig_cols, i % fig_cols]
            sns.histplot(data=data, ax=subplot, bins=bins)
            self._set_histogram_ax_properties(subplot, max_pos, group, i, fig_cols)
        plt.savefig(
            f"{self.project_dir}/output/enriched_{self.k}-mers_coverage_histograms.svg",
            format="svg",
            dpi=300,
            bbox_inches="tight",
        )

    def save_mappings_df(self):
        print(f"\nSaving {self.k}-mer mappings.")
        df_to_save = self.matrices["mappings"].copy()
        df_to_save.index.name = "id"
        df_to_save.to_csv(
            f"{self.project_dir}/output/enriched_{self.k}mers_mappings.csv"
        )
