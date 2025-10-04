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
        self.kmer_index = {}
        self.matrices = {}

        try:
            self.matrices["enriched"] = pd.read_csv(
                f"{project_dir}/output/enriched_{self.k}mers_stats.csv", index_col=0
            )
        except:
            raise FileNotFoundError(
                f"No enriched {self.k}-mers table found in {project_dir}output/ ! Please run PORT-EK find_enriched first!"
            )
        self.mapping_df = pd.DataFrame(
            [["", 0]],
            columns=["reference_sequence_position", "number_of_mismatches"],
            index=self.matrices["enriched"].index,
        )

    def index_ref_seq(self, max_mismatches: int, verbose: bool = False) -> None:
        if self._check_index(max_mismatches) == False:
            self._generate_and_write_index(max_mismatches, verbose)
        else:
            self._load_index(max_mismatches, verbose)

    def _check_index(self, max_mismatches: int) -> bool:
        return pathlib.Path(
            f"{self.project_dir}/temp/ref_index/{self.ref_seq_name}_index_{self.k}_{max_mismatches}.pkl"
        ).exists()

    def _generate_and_write_index(self, max_mismatches, verbose):
        if verbose == True:
            print(
                f"No {self.k}-mer index with maximum distance {max_mismatches} exists for {self.ref_seq_name}, building index."
            )
        kmer_index = self._generate_index(max_mismatches, self.ref_seq)
        os.makedirs(f"{self.project_dir}/temp/ref_index", exist_ok=True)
        with open(
            f"{self.project_dir}/temp/ref_index/{self.ref_seq_name}_index_{self.k}_{max_mismatches}.pkl",
            mode="wb",
        ) as out_file:
            pickle.dump(kmer_index, out_file)

        if verbose == True:
            print(
                f"Finished building {self.k}-mer index with maximum distance {max_mismatches} exists for {self.ref_seq_name}"
            )
            for d in range(max_mismatches + 1):
                print(f"Extracted {len(kmer_index[d])} k-mers with distance {d}.")
        self.kmer_index = kmer_index

    def _generate_index(self, max_mismatches: int, ref_seq: str) -> dict:
        kmer_index = {dist: {} for dist in range(max_mismatches + 1)}
        for kmer_start_position in range(0, len(ref_seq) - self.k + 1):
            kmer = ref_seq[kmer_start_position : kmer_start_position + self.k]
            self._save_kmer_to_index(kmer_index, kmer_start_position, 0, kmer)

            for n_mismatches in range(1, max_mismatches + 1):
                sub_kmers = self._generate_sub_kmers(n_mismatches, kmer)
                indel_kmers = self._generate_indel_kmers(
                    n_mismatches, kmer, kmer_start_position, ref_seq
                )
                ambi_kmers = sub_kmers.union(indel_kmers)
                for ambi_kmer in ambi_kmers:
                    self._save_kmer_to_index(
                        kmer_index, kmer_start_position, n_mismatches, ambi_kmer
                    )

        return kmer_index

    def _save_kmer_to_index(self, kmer_index, kmer_start_position, n_mismatches, kmer):
        if kmer in kmer_index[n_mismatches].keys():
            kmer_index[n_mismatches][kmer].append(kmer_start_position + 1)
        else:
            kmer_index[n_mismatches][kmer] = [kmer_start_position + 1]

    def _generate_sub_kmers(
        self,
        n_mismatches: int,
        kmer: str,
    ) -> set:

        sub_kmers = set()
        position_combinations = list(
            itertools.combinations(range(self.k), n_mismatches)
        )
        for positions in position_combinations:
            ambi_kmer = [nuc for nuc in kmer]
            for pos in positions:
                ambi_kmer[pos] = "N"
            sub_kmers.add("".join(ambi_kmer))

        return sub_kmers

    def _generate_indel_kmers(
        self,
        n_mismatches: int,
        kmer: str,
        kmer_start_position: int | None,
        ref_seq: str | None,
        generate_deletions: bool = True,
    ) -> set:

        indel_kmers = set()
        for position in range(1, self.k):
            in_kmer = [nuc for nuc in kmer]
            in_kmer.insert(position, "N" * n_mismatches)
            indel_kmers.add("".join(in_kmer))
            if generate_deletions == True:
                if position + n_mismatches >= self.k:
                    continue
                if kmer_start_position + self.k + n_mismatches > len(ref_seq):
                    continue
                del_kmer = [nuc for nuc in kmer]
                del_kmer = (
                    del_kmer[:position]
                    + del_kmer[position + n_mismatches :]
                    + [
                        nuc
                        for nuc in ref_seq[
                            kmer_start_position
                            + self.k : kmer_start_position
                            + self.k
                            + n_mismatches
                        ]
                    ]
                )
                indel_kmers.add("".join(del_kmer))

        return indel_kmers

    def _load_index(self, max_del_distance, verbose):
        if verbose == True:
            print(
                f"{self.k}-mer index with maximum distance {max_del_distance} for {self.ref_seq_name} already exists."
            )
        with open(
            f"{self.project_dir}/temp/ref_index/{self.ref_seq_name}_index_{self.k}_{max_del_distance}.pkl",
            mode="rb",
        ) as in_file:
            kmer_index = pickle.load(in_file)
        self.kmer_index = kmer_index

    def run_mapping(self, max_n_mismatch: int, verbose: bool = False) -> None:

        for kmer in self.mapping_df.index:
            mapping_dict = {n: set() for n in range(max_n_mismatch + 1)}
            for n_mismatch in range(max_n_mismatch + 1):
                self._map_kmer_to_index(kmer, n_mismatch, mapping_dict)
            for n_mismatch in range(max_n_mismatch + 1):
                if len(mapping_dict[n_mismatch]) > 0:
                    self.mapping_df.at[kmer, "reference_sequence_position"] = ",".join(
                        [str(pos) for pos in sorted(mapping_dict[n_mismatch])]
                    )
                    self.mapping_df.at[kmer, "number_of_mismatches"] = n_mismatch
                    break

    def _map_kmer_to_index(
        self, kmer: str, n_mismatches: int, mapping_dict: dict
    ) -> None:
        sub_kmers = self._generate_sub_kmers(n_mismatches, kmer)
        in_kmers = self._generate_indel_kmers(
            n_mismatches, kmer, None, None, generate_deletions=False
        )
        ambi_kmers = sub_kmers.union(in_kmers)
        for ambi_kmer in ambi_kmers:
            if ambi_kmer in self.kmer_index[n_mismatches].keys():
                mapping_dict[n_mismatches] = mapping_dict[n_mismatches].union(
                    self.kmer_index[n_mismatches][ambi_kmer]
                )

    def save_mapping(self, max_n_mismatch: int) -> None:
        self.mapping_df.to_csv(
            f"{self.project_dir}/output/mapping_{self.k}mers_max{max_n_mismatch}mismatches.csv"
        )
