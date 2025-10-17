import numpy as np
import pandas as pd
import networkx as nx
import regex
import os
import yaml
import pathlib
import pickle
from matplotlib import pyplot, collections, colormaps, patches, colors
from datetime import datetime
from scipy import stats
from scipy.spatial import distance
from scipy.cluster import hierarchy
from Bio import SeqIO, SeqRecord, Seq


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


def encode_seq_as_bits(seq: str) -> list[int]:

    encoding = {
        "A": 0b0001,  # A = 1
        "T": 0b0010,  # T = 2
        "G": 0b0100,  # G = 4
        "C": 0b1000,  # C = 8
        "W": 0b0011,  # W = A | T = 3
        "R": 0b0101,  # R = A | G = 5
        "Y": 0b1010,  # Y = C | T = 10
        "S": 0b1100,  # S = G | C = 12
        "K": 0b0110,  # K = G | T = 6
        "M": 0b1001,  # M = A | C = 9
        "B": 0b1110,  # B = C | G | T = 14
        "D": 0b0111,  # D = A | G | T = 7
        "H": 0b1011,  # H = A | C | T = 11
        "V": 0b1101,  # V = A | C | G = 13
        "N": 0b1111,  # N = A | C | G | T = 15
    }
    return [encoding.get(nuc.upper(), 15) for nuc in seq]


def encode_kmer(kmer_seq: str) -> int:
    encoding_dict = {"A": "00", "C": "01", "G": "10", "T": "11"}
    kmer_bin_string = [encoding_dict[nuc] for nuc in kmer_seq]
    id = int("".join(kmer_bin_string), base=2)
    return id


def decode_kmer(id: int, k) -> str:
    decoding = {"00": "A", "01": "C", "10": "G", "11": "T"}
    kmer_bin_string = bin(id)[2:].rjust(2 * k, "0")
    kmer_bin_string = [
        kmer_bin_string[i : i + 2] for i in range(0, len(kmer_bin_string), 2)
    ]
    kmer_seq = "".join([decoding[bits] for bits in kmer_bin_string])
    return kmer_seq


def filter_kmers(kmer_df: pd.DataFrame, freq_cols: list, cons_thr=0.01) -> pd.DataFrame:
    out_kmer_df = kmer_df[(kmer_df[freq_cols] > cons_thr).any(axis=1)]
    return out_kmer_df


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


# deprecated
def save_kmers_fasta(
    kmers: list, ids: list | None, name: str, directory: str, k: int
) -> None:
    out_fasta_list = []
    for i, kmer in enumerate(kmers):
        if ids == None:
            out_fasta_list.append(
                SeqRecord.SeqRecord(
                    seq=Seq.Seq(kmer),
                    id=f"{kmer}",
                    description="",
                )
            )
        else:
            out_fasta_list.append(
                SeqRecord.SeqRecord(
                    seq=Seq.Seq(kmer),
                    id=ids[i],
                    description="",
                )
            )
    SeqIO.write(
        out_fasta_list,
        f"{directory}/temp/{name}_{k}mers.fasta",
        format="fasta",
    )


# not implemented
def assemble_kmers_debruijn(kmers: list) -> nx.MultiDiGraph:
    de_bruijn = nx.MultiDiGraph()
    for kmer in kmers:
        kmerL = kmer[:-1]
        kmerR = kmer[1:]
        de_bruijn.add_node(kmerL)
        de_bruijn.add_node(kmerR)
        de_bruijn.add_edge(kmerL, kmerR, sequence=kmer)
    return de_bruijn


# not implemented
def assemble_contig(contig_graph: nx.MultiDiGraph) -> str:
    kmers = []
    if contig_graph.number_of_edges() == 1:
        path = nx.eulerian_path(contig_graph, keys=True)
        for edge in path:
            sequence = "".join([edge[0], edge[1][-1:]])
            kmers.append(nx.get_edge_attributes(contig_graph, "sequence")[edge])
    else:
        path = nx.eulerian_path(contig_graph, keys=True)
        first = True
        for edge in path:
            if first == True:
                sequence = "".join([edge[0][:1], edge[1]])
                first = False
            else:
                sequence = "".join([sequence, edge[1][-1:]])
            kmers.append(nx.get_edge_attributes(contig_graph, "sequence")[edge])
    return sequence, kmers


# not implemented
def cluster_kmers(
    matrix: pd.DataFrame, discard_singles: bool = False, max_d: float = 0
):
    distances = distance.pdist(matrix)
    linkage = hierarchy.linkage(distances)
    clustering = hierarchy.fcluster(linkage, max_d, criterion="distance")
    clustering = pd.Series(clustering, index=matrix.index, name="freq_cluster")
    freq_clusters = {cluster_no: [] for cluster_no in clustering.unique()}
    for cluster_no in freq_clusters.keys():
        freq_clusters[cluster_no] = clustering.loc[
            clustering == cluster_no
        ].index.to_list()
    cluster_graphs = {}
    for cluster, kmers in freq_clusters.items():
        cluster_graphs[cluster] = assemble_kmers_debruijn(kmers)
    contigs_kmer_dict = {}
    kmer_contig_dict = {}
    for cluster, graph in cluster_graphs.items():
        components = [
            graph.subgraph(component).copy()
            for component in nx.weakly_connected_components(graph)
        ]
        for component in components:
            if nx.is_directed_acyclic_graph(component):
                if nx.has_eulerian_path(component):
                    if discard_singles == True:
                        if component.number_of_edges() == 1:
                            continue
                    sequence, kmer_list = assemble_contig(component)
                    contigs_kmer_dict[sequence] = kmer_list
                    for kmer in kmer_list:
                        kmer_contig_dict[kmer] = sequence
                else:
                    print(f"Ambiguous assembly in cluster {cluster}")
            else:
                print(f"Ambiguous assembly in cluster {cluster}")
    return kmer_contig_dict, contigs_kmer_dict


# not_implemented
def make_ordinal(n):
    """
    Convert an integer into its ordinal representation::

        make_ordinal(0)   => '0th'
        make_ordinal(3)   => '3rd'
        make_ordinal(122) => '122nd'
        make_ordinal(213) => '213th'
    """
    n = int(n)
    if 11 <= (n % 100) <= 13:
        suffix = "th"
    else:
        suffix = ["th", "st", "nd", "rd", "th"][min(n % 10, 4)]
    return str(n) + suffix


# deprecated
def assemble_kmers(
    kmer_list: list, how: str = "pos", kmer_df: pd.DataFrame = None
) -> nx.DiGraph:
    edge_tuples = set()
    if how == "seq":
        for i in range(1, len(kmer_list)):
            for j in range(i):
                k1 = kmer_list[i]
                k2 = kmer_list[j]
                if k1[1:] == k2[:-1]:
                    edge_tuples.add((k1, k2))
                elif k2[1:] == k1[:-1]:
                    edge_tuples.add((k2, k1))
    elif how == "pos":
        if type(kmer_df) != pd.DataFrame:
            raise ValueError(
                "If assemblying kmers by position, you must provide a dataframe with kemrs as index and 'ref_start' column"
            )
        else:
            kmer_pos = [kmer_df.loc[kmer, "ref_start"] for kmer in kmer_list]
            for i in range(1, len(kmer_list)):
                for j in range(i):
                    k1 = kmer_list[i]
                    k2 = kmer_list[j]
                    if kmer_pos[i] == kmer_pos[j] - 1 and k1[1:] == k2[:-1]:
                        edge_tuples.add((k1, k2))
                    elif kmer_pos[j] == kmer_pos[i] - 1 and k2[1:] == k1[:-1]:
                        edge_tuples.add((k2, k1))

    kmer_graph = nx.from_edgelist(edge_tuples, create_using=nx.DiGraph)
    return kmer_graph


# not implemented
def plot_segments(segment_df, ref_seq, colormap=colormaps["coolwarm"]):
    plot_df = segment_df[["ref_start", "ref_end", "group"]].sort_values("ref_start")
    groups_values = ["ref"] + np.unique(plot_df["group"]).tolist()
    colors = [colormap(i) for i in np.linspace(0, 1, len(groups_values))]
    cdict = {groups_values[i]: colors[i] for i in range(len(groups_values))}
    ends_list = list(zip(plot_df["ref_start"], plot_df["ref_end"]))
    ys = []
    y = 1
    end_dict = {1: 0}
    while len(ends_list) > 1:
        ends = ends_list.pop(0)
        free_ends = {endy: endx for endy, endx in end_dict.items() if endx <= ends[0]}
        if len(free_ends.keys()) > 0:
            y = min(free_ends.keys())
        else:
            y = max(end_dict.keys()) + 1
        ys.append(y)
        end_dict[y] = ends[1]

    segments = [((1, 0), (len(ref_seq), 0))] + list(
        zip(zip(plot_df["ref_start"], ys), zip(plot_df["ref_end"], ys))
    )
    segment_colors = [cdict["ref"]] + [cdict[group] for group in plot_df["group"]]
    lines = collections.LineCollection(segments, colors=segment_colors, linewidths=2)
    fig, ax = pyplot.subplots()
    ax.set_xlim(0, len(ref_seq) + 1)
    ax.set_ylim(-1, max(ys) + 1)
    ax.add_collection(lines)
    legends = [patches.Patch(color=cdict[group], label=group) for group in cdict.keys()]
    pyplot.legend(handles=legends).set_loc("right")
    pyplot.show()


# not implemented
def _draw_genome_overlay_plot(
    segment_coords: list,
    segment_colors: list,
    ref_seq: str,
    title: str = None,
    colormap: colors.LinearSegmentedColormap = colormaps["coolwarm"],
    save_path: str = None,
    save_format: str = "svg",
):
    groups_values = ["ref"] + sorted(list(set(segment_colors)))
    colors = [colormap(i) for i in np.linspace(0, 1, len(groups_values))]
    cdict = {groups_values[i]: colors[i] for i in range(len(groups_values))}
    segment_colors = [cdict["ref"]] + [cdict[group] for group in segment_colors]
    segment_coords = [((1, 0), (len(ref_seq), 0))] + segment_coords
    lines = collections.LineCollection(
        segment_coords, colors=segment_colors, linewidths=2
    )
    fig, ax = pyplot.subplots()
    pyplot.tight_layout()
    ax.add_collection(lines)
    ax.autoscale()
    pyplot.tick_params(labelleft=False, left=False)
    ax.set_xlabel("Reference genome position")
    ax.set_ylabel("Genomes")
    ax.set_title(title)
    legends = [patches.Patch(color=cdict[group], label=group) for group in cdict.keys()]
    legends.reverse()
    pyplot.legend(handles=legends, bbox_to_anchor=(1.0, 0.5), loc="center left")
    if save_path != None:
        pyplot.savefig(save_path, format=save_format, dpi=600, bbox_inches="tight")
    pyplot.show()


# not implemented
def plot_kmers_by_genome(
    kmers_and_genomes: list,
    kmer_matrix: pd.DataFrame,
    group_sample_id: dict,
    ref_seq: str,
    title: str = None,
    colormap: colors.LinearSegmentedColormap = colormaps["coolwarm"],
    save_path: str = None,
    save_format: str = "svg",
):
    segment_coords = []
    segment_colors = []
    y = 1
    for kmers, samples in kmers_and_genomes:
        samples_to_plot = samples
        kmers_to_plot = kmers
        for sample_name in samples_to_plot:
            sample_group = [
                key for key, value in group_sample_id.items() if sample_name in value
            ][0]
            for kmer in kmers_to_plot:
                if kmer_matrix.loc[kmer, sample_name] > 0:
                    temp_segments = kmer_matrix.loc[kmer, "ref_pos"]
                    temp_segments = [((x1, y), (x2, y)) for x1, x2 in temp_segments]
                    segment_coords.extend(temp_segments)
                    segment_colors.extend(
                        [sample_group for _ in range(len(temp_segments))]
                    )
            y += 1

    _draw_genome_overlay_plot(
        segment_coords, segment_colors, ref_seq, title, colormap, save_path, save_format
    )


# not implemented
def assign_gene_from_interval(ref_pos: list, gene_dict: dict) -> str:
    genes = []
    for start, end in ref_pos:
        for gene, gene_ranges in gene_dict.items():
            for gene_range in gene_ranges:
                if (
                    len(
                        [
                            pos
                            for pos in range(start, end + 1)
                            if pos in list(range(gene_range[0], gene_range[1] + 1))
                        ]
                    )
                    > 0
                ):
                    genes.append(gene)


# not implemented
def assign_gene_from_position(ref_pos: int, gene_dict: dict) -> str:
    genes = []
    for gene, gene_ranges in gene_dict.items():
        for gene_range in gene_ranges:
            if gene_range[0] < ref_pos < gene_range[1]:
                genes.append(gene)
    return ", ".join(genes)


# deprecated
def calc_kmer_pvalue(kmer: str, first_group, sec_group, matrix: pd.DataFrame):
    first_obs = matrix.loc[kmer, first_group]
    sec_obs = matrix.loc[kmer, sec_group]
    test_result = stats.mannwhitneyu(first_obs, sec_obs)
    return test_result.pvalue


# Deprecated
def build_similarity_graph_two_list(
    name, query_list, target_list, mismatch_treshold: int
) -> nx.Graph:

    tot = len(query_list)
    similarity_edges = set()
    i = 1
    t0 = datetime.now()

    for query_kmer in query_list:
        similarity_edges.add((query_kmer, query_kmer))
        for target_kmer in target_list:
            if (
                sum(c1 != c2 for c1, c2 in zip(query_kmer, target_kmer))
                <= mismatch_treshold
            ):
                similarity_edges.add((query_kmer, target_kmer))
        print(
            f"Done {i} out of {tot} kmers. Time per kmer: {datetime.now()-t0}",
            end="\r",
            flush=True,
        )
        t0 = datetime.now()
        i += 1

    silmilarity_graph = nx.Graph(similarity_edges)
    return name, silmilarity_graph


# Deprecated
def calc_agg_freq(kmer_list, sample_list, source_df):
    pos_samples = set()
    for sample in sample_list:
        for kmer in kmer_list:
            if source_df.loc[kmer, sample] > 0:
                pos_samples.add(sample)
                break
    agg_freq = len(pos_samples) / len(sample_list)
    return agg_freq


# Deprecated
def map_kmers_find_mutations(kmer, ref_seq_str, pos_matrix, n=2, l=1000, find_wt=False):
    for i in range(n + 1):
        m = regex.findall(f"({kmer}){{s<={i}}}", ref_seq_str)
        if len(m) != 0:
            break
    if len(m) != 0:
        idx = [ref_seq_str.index(seq) + 1 for seq in m]
        mean_pos = pos_matrix.loc[kmer].mean()
        offset = list(abs(idx - mean_pos))

        if min(offset) <= l:
            best_align = offset.index(min(offset))
            best_match = m[best_align]
            start = idx[best_align]
            matchlist = [c1 == c2 for c1, c2 in zip(kmer, best_match)]
            n_match = sum(matchlist)
            alignment = {"seq": kmer, "start": start, "n_match": n_match}

            mutations = {"id": [], "ref_nt": [], "pos": [], "mut_nt": [], "kmer": []}
            for i, is_match in enumerate(matchlist):
                if find_wt == False:
                    if is_match == False:
                        mutations["ref_nt"].append(best_match[i])
                        mutations["pos"].append(start + i)
                        mutations["mut_nt"].append(kmer[i])
                        mutations["kmer"].append(kmer)
                        mutations["id"].append(f"{best_match[i]}{start+i}{kmer[i]}")

                if find_wt == True:
                    if is_match == True:
                        mutations["ref_nt"].append(best_match[i])
                        mutations["pos"].append(start + i)
                        mutations["mut_nt"].append(kmer[i])
                        mutations["kmer"].append(kmer)
                        mutations["id"].append(f"{best_match[i]}{start+i}{kmer[i]}")

        else:
            alignment = None
            mutations = {"id": [], "ref_nt": [], "pos": [], "mut_nt": [], "kmer": []}

    else:
        alignment = None
        mutations = {"id": [], "ref_nt": [], "pos": [], "mut_nt": [], "kmer": []}
    return alignment, mutations
