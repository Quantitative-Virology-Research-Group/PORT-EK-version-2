import concurrent.futures
import itertools
import pathlib
import pickle

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn import decomposition

import portek
from portek.portek_utils import BasePipeline


class EnrichedKmersPipeline(BasePipeline):
    """
    EnrichedKmersPipeline:
    """

    def __init__(self, project_dir: str, k: int) -> None:
        super().__init__(project_dir, k)
        super().load_sample_list()
        super().load_kmer_set(k)
        self.freq_cols = [f"{group}_freq" for group in self.sample_groups]
        self.avg_cols = [f"{group}_avg" for group in self.sample_groups]
        self.c_cols = [f"{group}_c" for group in self.sample_groups]
        self.f_cols = [f"{group}_f" for group in self.sample_groups]
        self.enriched_groups = None
        self.err_cols = None
        self.p_cols = None
        self.matrices = {}
        print(
            f"\nImported {len(self.kmer_set)} kmers and {len(self.sample_list)} samples."
        )

    def _filter_by_entropy_and_freq(
        self, matrix: pd.DataFrame, max_mem: float = 2.0, min_freq: float = 0.1
    ) -> pd.DataFrame:
        tot_samples = len(self.sample_list)
        matrix["F"] = matrix[self.c_cols].sum(axis=1) / tot_samples

        with np.errstate(divide="ignore"):
            matrix["H"] = np.where(
                matrix["F"] == 1,
                0,
                -(
                    matrix["F"] * np.log2(matrix["F"])
                    + (1 - matrix["F"]) * np.log2(1 - matrix["F"])
                ),
            )

        min_F = 2 / tot_samples
        min_H = -(min_F * np.log2(min_F) + (1 - min_F) * np.log2(1 - min_F))
        common_kmer_matrix = matrix.loc[
            (matrix["H"] >= min_H) | (matrix["F"] >= (1 - min_F))
        ]
        non_singles = len(common_kmer_matrix)

        print(f"{non_singles} {self.k}-mers passed the entropy filter.")

        if non_singles * len(self.sample_list) > max_mem * (2**30):
            common_kmer_matrix = common_kmer_matrix[
                ((common_kmer_matrix[self.freq_cols] > min_freq).sum(axis=1)) > 0
            ]
            print(
                f"The resulting count matrix would take over {max_mem} GB of memory. Removing additional {non_singles-len(common_kmer_matrix)} rare {self.k}-mers."
            )
            print(f"{len(common_kmer_matrix)} {self.k}-mers remaining.")

        return common_kmer_matrix

    def get_basic_kmer_stats(self, max_mem: int = 2, min_freq: float = 0.1) -> None:
        all_kmer_matrix = pd.DataFrame(
            0.0,
            index=self.kmer_set,
            columns=self.freq_cols + self.avg_cols,
            dtype=np.float64,
        ).sort_index()

        group_len_dict = {
            f"{group}": len(self.sample_group_dict[group])
            for group in self.sample_groups
        }

        in_path = pathlib.Path(f"{self.project_dir}/input/indices/").glob(
            f"{self.k}mer_*_avg_dict.pkl"
        )
        for filename in in_path:
            with open(filename, mode="rb") as in_file:
                temp_dict = pickle.load(in_file)
            column_name = "_".join(filename.stem.split("_")[1:-1])
            all_kmer_matrix[column_name] = temp_dict

        in_path = pathlib.Path(f"{self.project_dir}/input/indices/").glob(
            f"{self.k}mer_*_freq_dict.pkl"
        )
        for filename in in_path:
            with open(filename, mode="rb") as in_file:
                temp_dict = pickle.load(in_file)
            column_name = "_".join(filename.stem.split("_")[1:-1])
            all_kmer_matrix[column_name] = temp_dict

        all_kmer_matrix = all_kmer_matrix.fillna(0.0)

        for c_col, freq_col, group in zip(
            self.c_cols, self.freq_cols, self.sample_groups
        ):
            all_kmer_matrix[c_col] = round(
                all_kmer_matrix[freq_col] * group_len_dict[group], 0
            )

        common_kmer_matrix = self._filter_by_entropy_and_freq(
            all_kmer_matrix, max_mem, min_freq
        )
        self.matrices["common"] = common_kmer_matrix

    def _compare_group_pair(
        self,
        matrix_type: str,
        group1: str,
        group2: str,
        verbose: bool = False,
        false_discovery_control: bool = False,
        chunk_size: int = 2000,
    ) -> tuple[str, str, pd.Series, pd.Series, pd.Series]:
        try:
            if verbose:
                print(f"Calculating differences for groups {group1} and {group2}.")
            kmer_index = self.matrices[matrix_type].index
            avg_counts_i = self.matrices[matrix_type][f"{group1}_avg"]
            avg_counts_j = self.matrices[matrix_type][f"{group2}_avg"]
            errors = avg_counts_i - avg_counts_j
            group1_samples = self.sample_group_dict[group1]
            group2_samples = self.sample_group_dict[group2]
            g1_data = self.matrices[matrix_type][group1_samples].to_numpy()
            g2_data = self.matrices[matrix_type][group2_samples].to_numpy()
            n_kmers = g1_data.shape[0]
            p_vals = np.empty(n_kmers, dtype=np.float64)
            for start in range(0, n_kmers, chunk_size):
                end = min(start + chunk_size, n_kmers)
                p_vals[start:end] = stats.mannwhitneyu(
                    g1_data[start:end], g2_data[start:end], axis=1
                ).pvalue
            p_values = pd.Series(p_vals, index=kmer_index)
            if false_discovery_control:
                p_values = pd.Series(
                    stats.false_discovery_control(p_values, method="by"),
                    index=kmer_index,
                )
                if verbose:
                    print(
                        f"Applied Benjamini-Yekutieli false discovery control to p-values for groups {group1} and {group2}."
                    )
            with np.errstate(divide="ignore"):
                log_p_values = pd.Series(-np.log10(p_values), index=kmer_index)
            if verbose:
                print(f"Done calculating differences for groups {group1} and {group2}.")
            return group1, group2, errors, p_values, log_p_values
        except Exception as e:
            raise RuntimeError(
                f"Error comparing groups {group1} and {group2}: {e}"
            ) from e

    def calc_kmer_stats(
        self,
        matrix_type: str,
        n_jobs: int = 4,
        verbose: bool = False,
        false_discovery_control: bool = False,
        mwu_chunk_size: int = 2000,
    ) -> None:
        print(f"\nGetting {matrix_type} {self.k}-mer counts.")
        kmer_index = self.matrices[matrix_type].index
        filenames = sorted(
            pathlib.Path(f"{self.project_dir}/input/indices/{self.k}mers").glob(
                "*_count.pkl"
            )
        )
        tot_files = len(filenames)

        # Pre-allocate the exact final shape — no intermediate union of all k-mers
        count_array = np.zeros((len(kmer_index), len(self.sample_list)), dtype=np.uint8)
        kmer_to_row = {kmer: i for i, kmer in enumerate(kmer_index)}
        sample_to_col = {name: j for j, name in enumerate(self.sample_list)}

        def _fill_column(filename):
            with open(filename, mode="rb") as f:
                d = pickle.load(f)
            name = "_".join(filename.stem.split("_")[:-1])
            col = sample_to_col.get(name)
            if col is None:
                return
            rows, vals = [], []
            for k, v in d.items():
                row = kmer_to_row.get(k)
                if row is not None:
                    rows.append(row)
                    vals.append(v)
            if rows:
                count_array[rows, col] = vals

        loaded = 0
        with concurrent.futures.ThreadPoolExecutor(max_workers=n_jobs) as executor:
            for _ in executor.map(_fill_column, filenames):
                loaded += 1
                if verbose:
                    print(
                        f"Loaded {self.k}-mers from {loaded} of {tot_files} samples.",
                        end="\r",
                        flush=True,
                    )

        count_df = pd.DataFrame(count_array, index=kmer_index, columns=self.sample_list)
        self.matrices[matrix_type] = pd.concat(
            [count_df, self.matrices[matrix_type]], axis=1
        )
        print(f"\nIdentifying enriched {self.k}-mers.")

        if self.mode == "ava":
            err_cols = []
            p_cols = []
            group_pairs = []
            for j in range(1, len(self.sample_groups)):
                for i in range(j):
                    group_pairs.append(
                        (
                            matrix_type,
                            self.sample_groups[i],
                            self.sample_groups[j],
                            verbose,
                            false_discovery_control,
                            mwu_chunk_size,
                        )
                    )

            results = [self._compare_group_pair(*args) for args in group_pairs]

            new_cols = {}
            for result in results:
                err_name = f"{result[0]}-{result[1]}_err"
                p_name = f"{result[0]}-{result[1]}_p-value"
                err_cols.append(err_name)
                p_cols.append(p_name)
                new_cols[err_name] = result[2]
                new_cols[p_name] = result[3]
                new_cols[f"-log10_{p_name}"] = result[4]

            err_arrays = np.stack([new_cols[c].to_numpy() for c in err_cols], axis=1)
            new_cols["RMSE"] = pd.Series(
                np.sqrt(np.mean(err_arrays**2, axis=1)), index=kmer_index
            )
            self.matrices[matrix_type] = pd.concat(
                [self.matrices[matrix_type], pd.DataFrame(new_cols, index=kmer_index)],
                axis=1,
            )
            self.matrices[matrix_type] = self.matrices[matrix_type].sort_values(
                "RMSE", ascending=False
            )
            group_series = self.matrices[matrix_type].apply(
                portek.assign_kmer_group_ava,
                p_cols=p_cols,
                avg_cols=self.avg_cols,
                freq_cols=self.freq_cols,
                err_cols=err_cols,
                axis=1,
            )
            excl_series = self.matrices[matrix_type].apply(
                portek.check_exclusivity, avg_cols=self.avg_cols, axis=1
            )
            self.matrices[matrix_type] = pd.concat(
                [
                    self.matrices[matrix_type],
                    pd.DataFrame({"group": group_series, "exclusivity": excl_series}),
                ],
                axis=1,
            )

        elif self.mode == "ovr":
            err_cols = []
            p_cols = []
            group_pairs = []
            for group in self.control_groups:
                group_pairs.append(
                    (
                        matrix_type,
                        self.goi,
                        group,
                        verbose,
                        false_discovery_control,
                        mwu_chunk_size,
                    )
                )

            results = [self._compare_group_pair(*args) for args in group_pairs]

            new_cols = {}
            for result in results:
                err_name = f"{result[0]}-{result[1]}_err"
                p_name = f"{result[0]}-{result[1]}_p-value"
                err_cols.append(err_name)
                p_cols.append(p_name)
                new_cols[err_name] = result[2]
                new_cols[p_name] = result[3]
                new_cols[f"-log10_{p_name}"] = result[4]

            err_arrays = np.stack([new_cols[c].to_numpy() for c in err_cols], axis=1)
            new_cols["RMSE"] = pd.Series(
                np.sqrt(np.mean(err_arrays**2, axis=1)), index=kmer_index
            )
            self.matrices[matrix_type] = pd.concat(
                [self.matrices[matrix_type], pd.DataFrame(new_cols, index=kmer_index)],
                axis=1,
            )
            self.matrices[matrix_type] = self.matrices[matrix_type].sort_values(
                "RMSE", ascending=False
            )
            group_series = self.matrices[matrix_type].apply(
                portek.assign_kmer_group_ovr,
                goi=self.goi,
                p_cols=p_cols,
                err_cols=err_cols,
                freq_cols=self.freq_cols,
                axis=1,
            )
            excl_series = self.matrices[matrix_type].apply(
                portek.check_exclusivity, avg_cols=self.avg_cols, axis=1
            )
            self.matrices[matrix_type] = pd.concat(
                [
                    self.matrices[matrix_type],
                    pd.DataFrame({"group": group_series, "exclusivity": excl_series}),
                ],
                axis=1,
            )

        self.enriched_groups = [
            name.split("_")[0]
            for name in self.matrices[matrix_type]["group"].value_counts().index
            if "enriched" in name
        ]
        self.err_cols = err_cols
        self.p_cols = p_cols

    def plot_volcanos(self, matrix_type):
        print(f"\nPlotting and saving volcano plots of enriched {self.k}-mers.")
        for i in range(len(self.err_cols)):
            err = self.err_cols[i]
            group1 = err.split("_")[0].split("-")[0]
            group2 = err.split("_")[0].split("-")[1]
            logp = f"-log10_{self.p_cols[i]}"
            fig, ax = plt.subplots()
            plt.axhline(y=-np.log10(0.01), color="black")
            plt.axvline(x=0.1, color="black", linestyle="--")
            plt.axvline(x=-0.1, color="black", linestyle="--")
            ax.autoscale()
            ax.set_title(f"{group1} vs {group2} volcano plot")
            ax.set_xlabel("Average kmer count change")
            ax.set_ylabel("-log10 of p-value")
            fig.tight_layout()
            sns.scatterplot(
                data=self.matrices[matrix_type],
                x=err,
                y=logp,
                s=10,
                linewidth=0,
                hue="group",
                alpha=0.5,
            )
            plt.savefig(
                f"{self.project_dir}/output/{err}_{matrix_type}_{self.k}mers_volcano.svg",
                dpi=600,
                format="svg",
                bbox_inches="tight",
            )

    def get_enriched_kmers(self):
        if "rare_similar" not in self.matrices.keys():
            self.matrices["enriched"] = self.matrices["common"].loc[
                (
                    (self.matrices["common"]["group"] != "not_significant")
                    & (self.matrices["common"]["group"] != "group_dependent")
                    & (self.matrices["common"]["RMSE"] > 0.1)
                )
                | (self.matrices["common"]["group"] == "conserved")
            ]
        else:
            self.matrices["enriched"] = pd.concat(
                [
                    self.matrices["common"].loc[
                        (
                            (self.matrices["common"]["group"] != "not_significant")
                            & (self.matrices["common"]["group"] != "group_dependent")
                            & (self.matrices["common"]["RMSE"] > 0.1)
                        )
                        | (self.matrices["common"]["group"] == "conserved")
                    ],
                    self.matrices["rare_similar"].loc[
                        (
                            (
                                self.matrices["rare_similar"]["group"]
                                != "not_significant"
                            )
                            & (
                                self.matrices["rare_similar"]["group"]
                                != "group_dependent"
                            )
                            & (self.matrices["rare_similar"]["RMSE"] > 0.1)
                        )
                        | (self.matrices["common"]["group"] == "conserved")
                    ],
                ]
            )
        if (
            len(
                self.matrices["enriched"][
                    self.matrices["enriched"]["group"] != "conserved"
                ]
            )
            == 0
        ):
            print(f"\nPORT-EK has found no enriched {self.k}-mers!")
            return False
        else:
            print("\nPORT-EK has found:")
            group_numbers = self.matrices["enriched"]["group"].value_counts()
            for group in group_numbers.index:
                print(f"{group_numbers.loc[group]} {group} {self.k}-mers")
            return True

    def plot_PCA(self):
        print(f"\nPlotting and saving PCA plot of enriched {self.k}-mers.")
        pca = decomposition.PCA(2)
        X_PCA = pca.fit_transform(self.matrices["counts"].drop("sample_group", axis=1))
        y_names = self.matrices["counts"]["sample_group"]
        fig, ax = plt.subplots()
        fig.tight_layout()
        ax.set_title("PCA of k-mer counts")
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
        sns.scatterplot(x=X_PCA[:, 0], y=X_PCA[:, 1], hue=y_names, s=20, linewidth=0)
        plt.savefig(
            f"{self.project_dir}/output/{self.k}mer_PCA.svg",
            dpi=600,
            format="svg",
            bbox_inches="tight",
        )

    def save_counts_for_classifier(self):
        print(f"\nSaving {self.k}-mer counts for classification.")
        counts_for_classifier = (
            self.matrices["enriched"]
            .loc[self.matrices["enriched"]["group"] != "conserved", self.sample_list]
            .T
        )
        counts_for_classifier.columns = counts_for_classifier.columns.map(
            lambda id: portek.decode_kmer(id, self.k)
        )
        counts_for_classifier["sample_group"] = counts_for_classifier.index.map(
            lambda name: name.split("_")[0]
        )
        self.matrices["counts"] = counts_for_classifier
        counts_for_classifier.to_csv(
            f"{self.project_dir}/output/{self.k}mer_counts_for_classifier.csv",
            index_label="sample_name",
        )

    def save_matrix(self, matrix_type: str, full: bool = False):
        print(f"\nSaving {matrix_type} {self.k}-mers matrix.")
        df_to_save = self.matrices[matrix_type].copy()
        df_to_save.index = df_to_save.index.map(
            lambda id: portek.decode_kmer(id, self.k)
        )
        if full == True:
            out_filename = f"{self.project_dir}/output/{matrix_type}_{self.k}mers.csv"
            df_to_save.to_csv(out_filename, index_label="kmer")
        else:
            out_filename = (
                f"{self.project_dir}/output/{matrix_type}_{self.k}mers_stats.csv"
            )
            export_cols = (
                self.avg_cols
                + list(itertools.chain(*zip(self.err_cols, self.p_cols)))
                + ["RMSE", "group", "exclusivity"]
            )
            df_to_save.loc[:, export_cols].to_csv(out_filename, index_label="kmer")
