#!/usr/bin/env python3
"""R1-Major3 / minor 3: per-group cohort description from LANL FASTA headers (subtype.country.year.name.accession).
usage: python s06_cohort_summary.py DATA_DIR DNA|RNA OUT_PREFIX
Reports n, exact duplicates (identical ungapped sequence), countries, years, accession-stem proxy for study/submission,
gap fraction and ambiguous-base fraction. Patient ID/study are NOT in the header: add them from the LANL metadata download
(Sequence Search interface -> 'Patient ID', 'Reference') and use the one-per-patient filter for deduplication."""
import sys, re, hashlib, pandas as pd
from Bio import SeqIO
d, typ, out = sys.argv[1:4]
files = {"A": "A.comp", "B": "B", "C": "C", "D": "D", "Rare": "Rare", "CRF01.AE": "CRF01.AE", "CRF02.AG": "CRF02.AG"}
rows = []
for g, f in files.items():
    for r in SeqIO.parse(f"{d}/hiv-db.{typ}.{f}.fasta", "fasta"):
        p = r.id.split("."); s = str(r.seq).upper(); u = re.sub("[^ACGT]", "", s)
        rows.append(dict(group=g, id=r.id, lanl_subtype=p[0], country=p[1], year=p[2], acc=p[-1], stem=re.sub(r"\d{2}$", "", p[-1]),
                         aln_len=len(s), nt_len=len(u), gap_frac=s.count("-") / len(s), amb_frac=sum(c not in "ACGT-" for c in s) / max(1, len(u)),
                         md5=hashlib.md5(u.encode()).hexdigest()))
df = pd.DataFrame(rows); df.to_csv(out + "_sequences.csv", index=False)
yr = pd.to_numeric(df.year, errors="coerce")
tab = df.assign(yr=yr).groupby("group", sort=False).agg(n=("id", "size"), unique_seq=("md5", "nunique"), countries=("country", "nunique"),
        year_min=("yr", "min"), year_max=("yr", "max"), submissions_proxy=("stem", "nunique"), nt_len=("nt_len", "median"),
        gap_frac=("gap_frac", "median"), amb_frac=("amb_frac", "median"))
tab["exact_dups"] = tab.n - tab.unique_seq; tab["top_country"] = df.groupby("group", sort=False).country.agg(lambda x: "; ".join(f"{k}:{v}" for k, v in x.value_counts().head(3).items()))
tab.round(3).to_csv(out + "_summary.csv"); print(tab.round(3).to_string())
