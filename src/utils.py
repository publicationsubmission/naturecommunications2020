#!/usr/bin/env python

import pandas as pd

from collections import defaultdict
from Bio import AlignIO


def read_alignment(fname: str, mapping: dict) -> dict:
    msa = AlignIO.read(fname, "fasta")

    seqs = {}

    for alignment in msa:
        seqs[alignment.id] = str(alignment.seq)

    return seqs


def tabulate_mutations(
    seqs: dict, mapping: dict, df: pd.DataFrame, output="output"
) -> pd.DataFrame:

    result = defaultdict(list)
    inverse_mapping = {v: k for k, v in mapping.items()}
    ref = list(filter(lambda x: not x.startswith("E_"), seqs.keys()))[0]
    pos = list(range(len(seqs[ref])))
    ref_bases = list(seqs[ref])
    df["Position"] = list(map(lambda x: x + 1, pos))
    df["Reference"] = ref_bases

    for sample in seqs.keys():
        if sample != ref:
            sample_bases = list(seqs[sample])
            df[mapping[sample]] = sample_bases

    df.fillna("", inplace=True)

    df.to_csv(f"{output}/{mapping[ref]}.csv", index=False)
    return df
